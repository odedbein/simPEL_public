function simPEL9_recognition(subj_num, subj_id, fMRI,practice, debug)
% function recognition(subj_num, subj_initials, fMRI, debug)
% Inputs:
%   subjectNumber: start first subject as 1, used to counterbalance button
%                       box response mappings (integer)
%   subj_initials:      subject initials, used to write output file names (string)
%   fMRI:               [0 / 1] looks for scanner pulse/trigger if fMRI,
%                       restricts KbCheck to button box
%   debug:              [0 / 1] if debug, not whole screen

%% boiler plate
if ~debug;
    ListenChar(2);
    HideCursor;
end

Priority(0);
Screen('Preference', 'SkipSyncTests', 1 );
description = Screen('Computer'); %#ok<NASGU>

%% paths
project_dir='.';
addpath(fullfile(project_dir,'scripts'))
stim_dir=fullfile(project_dir,'images');
data_dir=fullfile(project_dir,'data',sprintf('%d%s',subj_num,subj_id));


if ~exist(data_dir,'dir')
    correct_part='n';
    while correct_part=='n'
        prompt='data_dir does not exist, are you sure you''re running the correct section?[y/terminate]';
        correct_part=input(prompt,'s');
    end
end

%% file lists:
categ_names={'small','big'};
file_names={};
for c=1:numel(categ_names)
    
    %that's for the foils:
    dir_list=dir([stim_dir '/' categ_names{c}]);
    % skip stuff that isn't an image
    dir_list = dir_list(3:end);               % skip . & ..
    if (strcmp(dir_list(1).name,'.DS_Store')==1) % also skip .DS_Store
        dir_list = dir_list(2:end);
    end
    for image=1:numel(dir_list)
        file_names{image,c}=[stim_dir '/' categ_names{c} '/' dir_list(image).name];
    end
    
    %that's for the paired
    dir_list=dir([stim_dir '/paired/' categ_names{c}]);
    % skip stuff that isn't an image
    dir_list = dir_list(3:end);               % skip . & ..
    if (strcmp(dir_list(1).name,'.DS_Store')==1) % also skip .DS_Store
        dir_list = dir_list(2:end);
    end
    for image=1:numel(dir_list)
        file_names{image,c+numel(categ_names)}=[stim_dir '/paired/' categ_names{c} '/' dir_list(image).name];
    end
end

%% experiment parameters
% stimulus and sequence info loaded in from encoding task

run_delay = 12;
run_decay = 12;
% end
if debug
    stim_dur = 0.1;           %sec
    response_window = 0.1;    %sec
    run_delay = 1;
    run_decay = 1;
else
    stim_dur = 3;           %sec
    response_window = 5;    %sec
end

%% keyboard and key presses

device=-1; %allows use all keyboard platforms
KbName('UnifyKeyNames');

if fMRI == 1
    LEFT = KbName('4$');
    RIGHT = KbName('3#');
elseif fMRI == 0
    LEFT = KbName('j');
    MIDDLE = KbName('k');
    RIGHT = KbName('l');
end

% counterbalance response presses by subject
if (mod(subj_num,2) == 0)
    instructString{1} = 'Decide whether each object is old, similar or new';
    instructString{2} = 'If you think that the object is old, press ''J'' ';
    instructString{3} = 'If you think that the object is similar, press ''K'' ';
    instructString{4} = 'If you think that the object is new, press ''L'' ';
    instructString{5} = 'Press ''q'' to begin!';
    
    old=LEFT; similar=MIDDLE; new=RIGHT;
else
    instructString{1} = 'Decide whether each object is old, similar or new';
    instructString{2} = 'If you think that the object is old, press ''L'' ';
    instructString{3} = 'If you think that the object is similar, press ''K'' ';
    instructString{4} = 'If you think that the object is new, press ''J'' ';
    instructString{5} = 'Press ''q'' to begin!';
    
    old=RIGHT; similar=MIDDLE; new=LEFT;
end

if practice == 1
    instructString{end+1} = '(practice)';
end

% instructions to get the experimenter at the end
if fMRI==0
    progressString{1} = 'Great job! Please find the';
    progressString{2} = 'experimenter in room 871.';
elseif fMRI==1
    progressString{1} = 'Great job!';
end

%confidence buttons:
%confInst{1} = 'OLD SURE        OLD UNSURE'; 
%confInst{2} = 'NEW UNSURE        NEW SURE';

%% set screen information
backColor = 127;
textColor = 0;
textFont = 'Arial';
textSize = 18;
textSpacing = 25;
fixColor = 0;
imageSize=350; % assumed square
fixationSize = 4; % pixels
screens = Screen('Screens');

if fMRI
    % present stimuli in second screen
    windowIdx=max(screens);
    [Screen_X, Screen_Y]=Screen('WindowSize',windowIdx);
elseif debug
    % present stimuli on second screen, but smaller
    windowIdx=max(screens);
    Screen_X = 640;
    Screen_Y = 480;
else
    % behavioral: run in testing room
    windowIdx=min(screens);
    [Screen_X, Screen_Y]=Screen('WindowSize',windowIdx);
end
centerX = Screen_X/2;
centerY = Screen_Y/2;

% placeholder for images
imageRect = [0,0,imageSize,imageSize];

% position of images
centerRect = [centerX-imageSize/2,centerY-imageSize/2,centerX+imageSize/2,centerY+imageSize/2];

% position of fixation dot
fixDotRect = [centerX-fixationSize,centerY-fixationSize,centerX+fixationSize,centerY+fixationSize];

%% load images

% load priming trial sequences
load([data_dir '/rec_trial_sequences.mat'],'rec_trial_seq','rec_onsets');
n_recog_runs=size(rec_trial_seq,1);

% create main window
mainWindow = Screen(windowIdx,'OpenWindow',backColor,[0 0 Screen_X Screen_Y]);
Screen(mainWindow,'TextFont',textFont);
Screen(mainWindow,'TextSize',textSize);

%% practice
if practice==1
    
    % override number of encoding trials and stim folder
    n_recog_runs=1;
    prac_idx=[1 -2 13 14 -15 6 -7 8 9 10 11 12];
    condition_temp=[6 6 8 8 8 6 6 6 9 9 9 9];%chose 8/6 just to have same/similar, could have been 7/5 just as well
    condition=zeros(size(prac_idx));
    dir_list=dir([stim_dir '/practice']);
    
    % skip stuff that isn't an image
    dir_list = dir_list(3:end);               % skip . & ..
    if (strcmp(dir_list(1).name,'.DS_Store')==1) % also skip .DS_Store
        dir_list = dir_list(2:end);
    end
    
    prac_file_names={};
    % grab some objects
    rand_rec=randperm(length(prac_idx));
    for ff=1:length(prac_idx)
        f=rand_rec(ff);
        condition(ff)=condition_temp(f);
        prac_file_names{ff}=[stim_dir '/practice/' dir_list(abs(prac_idx(f))).name];
    end
    
end

%% grab the images
if practice==1
    prim_texture=zeros(size(prac_idx));
    for i=1:length(prac_idx) %the number of the stimuli
        
        %enc_idx is the order of the stimuli
        % match encoding index to the list of all indices for this subject
        %all_idx is an array of numbers - possibly, this is the array that sets
        %which image is "image number x" it is different in each subject
        %takes "image number 1" finds it in loc - say, number 40.
        %temp_idx=all_idx==i;
        
        % read in images in the order of encoding presentation
        temp_image = imread(prac_file_names{i}); %now in temp image, the first 1, there is the 40th image in the folder
        
        % return index of images that can be called by 'DrawTexture'
        prim_texture(i) = Screen('MakeTexture',mainWindow,temp_image);
        %makes an array of images ready to be displayed
        %put in where there is the number "1" in enc_idx, the 40th image in the folder.
        %so, all_idx and he order in which items appear in enc_idx set the
        %randomization
    end
    
    %these are unecessary in the practice, they are here only so that
    %the script won't break
    %preparing a stimuli size array
    priming_sequence=prac_idx(rand_rec);
else %not practice
    %grab currtn trial sequence, response,condition:
    priming_sequence=rec_trial_seq;
    condition=rem(abs(priming_sequence),10000);
    condition=floor(condition/1000);
    
    prim_texture=zeros(size(priming_sequence));
    for i=1:numel(priming_sequence)
        curr_image=priming_sequence(i);
        image_size=floor(abs(curr_image)/10000);
        image_num=rem(abs(curr_image),1000);
        % read in images in the order of encoding presentation
        if condition(i)==9 %take from the foils array
            temp_image = imread(file_names{image_num,image_size});
        else
            temp_image = imread(file_names{image_num,image_size+numel(categ_names)});
        end
        
        % return index of images that can be called by 'DrawTexture'
        prim_texture(priming_sequence==curr_image) = Screen('MakeTexture',mainWindow,temp_image);
        
    end
end

%set up some variables:
stim_onset=zeros(size(prim_texture));
response_type=zeros(size(prim_texture));
response_accuracy=zeros(size(prim_texture));
response_times=zeros(size(prim_texture));
response=zeros(size(prim_texture));

%% run each priming section:
for run=1:n_recog_runs

    %% trial onsets and data structures
    trial_onsets=rec_onsets(run,:);
    tpr=length(priming_sequence(run,:));
    
    if debug == 1
        iti= 0.1 * ones(size(trial_onsets));
        trial_lengths=stim_dur + iti;
        trial_onsets=run_delay + cumsum(trial_lengths) - trial_lengths;
    else
        if practice==1 %this is practice - do not wait the initial 12 secs
            run_delay = 1;
            run_decay = 1;
            rand_onsets=[1 1 3 5 1 1 3 5 1 1 3 5];
            rand_onsets=rand_onsets(randperm(length(rand_onsets)));
            trial_lengths=stim_dur + rand_onsets;
            trial_onsets=cumsum(trial_lengths,2) - (trial_lengths)+run_delay;
            %run_length=trial_onsets(end)+stim_dur+1; % end 1 sec after last trial
        end
    end
    
    
    run_length=trial_onsets(end) + stim_dur + run_decay;
    
    % set up output text file
    if practice==0
        data_file=fullfile(data_dir,['recognition_' num2str(subj_num) subj_id 'part_' num2str(run) '.txt']);
        if exist(data_file,'file')
            data_file=[data_file(1:(end-4)) '1.txt'];
        end
        fid = fopen(data_file, 'w');
        
        fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', ...
            'subj', 'init', 'run', 'trial', 'onset', 'image', 'condition', 'resp', 'resp_type', 'acc', 'rt');
    end
    % print header to screen
    fprintf('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', ...
        'subj', 'init', 'run', 'trial', 'onset', 'image', 'condition', 'resp', 'resp_type', 'acc', 'rt');

    % set up diagnostic plot
%     if fMRI==1 || debug==1
%         rd = response_diagnostics(response_window, tpr);
%     end

    %% show instructions
    for i=1:length(instructString)
        tempBounds = Screen('TextBounds',mainWindow,instructString{i});
        Screen('drawtext',mainWindow,instructString{i},round(centerX-tempBounds(3)/2),round(centerY-tempBounds(4)/5+textSpacing*(i-1)),textColor);
        clear tempBounds;
    end
    Screen('Flip',mainWindow);
    
    % wait for experimenter to advance with 'q' key
    FlushEvents('keyDown');
    while(1)
        temp = GetChar;
        if (temp == 'q')
            break;
        end
    end
    
    %clear screen
    Screen(mainWindow,'FillRect',backColor);
    Screen('Flip',mainWindow);
    runTime = GetSecs;
    
    % fixation until first trial
    Screen(mainWindow,'FillOval',fixColor,fixDotRect);
    Screen('Flip',mainWindow);
    
    %% trial loop
    
    for trial=1:tpr
        
        % wait for next trial onset
        while GetSecs<runTime+trial_onsets(trial); end
        
        % present images
        Screen('DrawTexture',mainWindow,prim_texture(run,trial),imageRect,centerRect);
%         tempBounds = Screen('TextBounds',mainWindow,confInst{1});
%         Screen('drawtext',mainWindow,confInst{1},round(Screen_X*.30-tempBounds(3)/2),round(Screen_Y*.80),textColor);
%         clear tempBounds;
%         tempBounds = Screen('TextBounds',mainWindow,confInst{2});
%         Screen('drawtext',mainWindow,confInst{2},round(Screen_X*.70-tempBounds(3)/2),round(Screen_Y*.80),textColor);
%         clear tempBounds;
        
        stim_start=Screen('Flip',mainWindow);
        stim_onset(run,trial)=stim_start-runTime;
        
        
        %  collect reponse
        removed=0;
        FlushEvents('keyDown');
        while GetSecs < stim_start + response_window
            WaitSecs(.0005);
            
            % take down image even if no response yet
            if  GetSecs >= (stim_start + stim_dur)
                %removed = 1;
                Screen(mainWindow,'FillOval',fixColor,fixDotRect);
%                 tempBounds = Screen('TextBounds',mainWindow,confInst{1});
%                 Screen('drawtext',mainWindow,confInst{1},round(Screen_X*.30-tempBounds(3)/2),round(Screen_Y*.80),textColor);
%                 clear tempBounds;
%                 tempBounds = Screen('TextBounds',mainWindow,confInst{2});
%                 Screen('drawtext',mainWindow,confInst{2},round(Screen_X*.70-tempBounds(3)/2),round(Screen_Y*.80),textColor);
%                 clear tempBounds;
                Screen('Flip', mainWindow);
            end
            
            % wait for response
            [keyIsDown, secs, keyCode] = KbCheck(device);
                if (keyIsDown)
                    %removed = 1;
                    % record response button, response time, and response
                    % confidence
                    response(run,trial) = find(keyCode,1);
                    response_times(run,trial) = secs - stim_start;
                    switch find(keyCode,1)
                        case old
                             response_type(run,trial)=3;
                        case similar
                             response_type(run,trial)=2;
                        case new
                             response_type(run,trial)=1;
                    end
                    % record accuracy
                    if condition(run,trial)==9 && (find(keyCode,1) == new)
                        response_accuracy(run,trial) = 1;
                    elseif (condition(run,trial)==5 || condition(run,trial)==6) && (find(keyCode,1) == old)
                        response_accuracy(run,trial) = 1;
                    elseif (condition(run,trial)==7 || condition(run,trial)==8) && (find(keyCode,1) == similar)
                        response_accuracy(run,trial) = 1;
                    else
                        response_accuracy(run,trial) = 0;
                    end
                    
                end
        end
        
        % fixation until next trial
        Screen(mainWindow,'FillOval',fixColor,fixDotRect);
%         tempBounds = Screen('TextBounds',mainWindow,confInst{1});
%         Screen('drawtext',mainWindow,confInst{1},round(Screen_X*.30-tempBounds(3)/2),round(Screen_Y*.80),textColor);
%         clear tempBounds;
%         tempBounds = Screen('TextBounds',mainWindow,confInst{2});
%         Screen('drawtext',mainWindow,confInst{2},round(Screen_X*.70-tempBounds(3)/2),round(Screen_Y*.80),textColor);
%         clear tempBounds;
        Screen('Flip',mainWindow);
        
        % print output to screen and text file
        fprintf('%d\t%s\t%d\t%d\t%.3f\t%d\t%d\t%d\t%d\t%d\t%.3f\n', ...
                subj_num, subj_id, run, trial, stim_onset(run,trial), priming_sequence(run,trial),condition(run,trial), ...
                response(run,trial),response_type(run,trial), response_accuracy(run,trial), response_times(run,trial));
            
        if practice==0
            fprintf(fid,'%d\t%s\t%d\t%d\t%.3f\t%d\t%d\t%d\t%d\t%d\t%.3f\n', ...
                        subj_num, subj_id, run, trial, stim_onset(run,trial), priming_sequence(run,trial),condition(run,trial), ...
                        response(run,trial), response_type(run,trial), response_accuracy(run,trial), response_times(run,trial));
            
        end
        
        % plot performance and rt ongoing
%         if fMRI==1 || debug==1
%             rd.log_trial(response_times(trial), correct_response(trial), response(trial));
%         end
    end
    
    % wait for scanner to finish
    if practice==0;
        while GetSecs < runTime + trial_onsets(end) + stim_dur + run_decay; end
        % close this run's text file
        fclose(fid);
    end
end

%% clean up and go home
if practice==0
    mat_file=fullfile(data_dir,'recognition.mat');
    if exist(mat_file,'file')
        mat_file=[mat_file(1:end-4) '1.mat'];
    end
    save(mat_file,'response','response_times','response_accuracy','response_accuracy','stim_onset','description','priming_sequence','condition','response_type');
    
    %% show instructions to get experimenter
    for i=1:length(progressString)
        tempBounds = Screen('TextBounds',mainWindow,progressString{i});
        Screen('drawtext',mainWindow,progressString{i},round(centerX-tempBounds(3)/2),round(centerY-tempBounds(4)/5+textSpacing*(i-1)),textColor);
        clear tempBounds;
    end
    Screen('Flip',mainWindow);
    
    % wait for experimenter to advance with 'q' key
    FlushEvents('keyDown');
    while(1)
        temp = GetChar;
        if (temp == 'q')
            break;
        end
    end

end

if ~debug
    ListenChar(0);
    ShowCursor;
end

% clear screen
sca