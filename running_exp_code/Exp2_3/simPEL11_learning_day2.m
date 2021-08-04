function simPEL11_learning_day2(subj_num, subj_id,varargin)
% function encoding(subj_num, subj_initials, fMRI, practice, debug, [startBlock])
% Inputs:
%   subj_num: start first subject as 1, used to counterbalance button
%                       box response mappings (integer)
%   subj_initials:      subject initials, used to write output file names (string)
%   part:               'day1_init'/'day2_init'/'rem'/'er_block'
%   fMRI:               [0 / 1] looks for scanner pulse/trigger if fMRI,
%                       restricts KbCheck to button box
%   practice:           [0 / 1] run through 10 practice trials
%   debug:              [0 / 1] if debug, not whole screen

%% boiler plate

fMRI=varargin{1};
practice=varargin{2};
debug=varargin{3};

if length(varargin)==4
    startBlock=varargin{4};
else
    startBlock=1;
end


if ~debug
    ListenChar(2);
    HideCursor;
end

Priority(0);
Screen('Preference', 'SkipSyncTests', 1 );
description = Screen('Computer'); %#ok<NASGU>

%make sure that epxerimenter inserted a correct name
%%%% what I really should do is to make it check there are available txt
%%%% files and shout if doesn't.
% correct_part='n';
% while correct_part=='n'
%     switch part
%         case 'day1_init'
%             prompt = 'are you sure that ''day1_init'' is the correct part?[y/terminate]';
%             correct_part=input(prompt,'s');
%         case 'day2_init'
%             prompt = 'are you sure that ''day2_init'' is the correct part?[y/terminate]';
%             correct_part=input(prompt,'s');
%         case 'rem'
%             prompt = 'are you sure that ''rem'' is the correct part?[y/terminate]';
%             correct_part=input(prompt,'s');
%         case 'er_block'
%             prompt = 'are you sure that ''er_blcok'' is the correct part?[y/terminate]';
%             correct_part=input(prompt,'s');
%         otherwise
%             display('wrong part, please terminate and restart the script');
%             while correct_part=='n'
%             end
%     end
% end
%% paths
project_dir='.';
addpath(fullfile(project_dir,'scripts'))
stim_dir=fullfile(project_dir,'images');
data_dir=fullfile(project_dir,'data',sprintf('%d%s',subj_num,subj_id));

if ~exist(data_dir,'dir')
    mkdir(data_dir);
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
%% trial timing: all units in sec
if debug
    stim_dur = 0.1;
    response_window = 0.1;
else
    stim_dur = 1.5;
    response_window = 1.9;
end
run_delay = 2;
run_decay = 2;

%% keyboard and key presses
device=-1; %allows use all keyboard platforms
KbName('UnifyKeyNames');

if fMRI == 1
    LEFT = KbName('4$');
    RIGHT = KbName('3#');
elseif fMRI == 0
    LEFT = KbName('j');
    RIGHT = KbName('k');
end

backtick=KbName('`~');

% counterbalance response presses by subject
%we ran two versions - one with item focused task (bigger/smaller than shoe box), one with associative
%focused task (bigger smaller than previous).
%Here is the item: instructString{1} = 'Decide whether each object is bigger or smaller than a shoe box in real life';
%that's the associative: instructString{1} = 'Decide whether each object is bigger or smaller than the previous object in real life';
if (mod(subj_num,2) == 0)
    instructString{1} = 'Decide whether each object is bigger or smaller than a shoe box in real life';
    instructString{2} = 'If bigger, press ''J'' (index finger)';
    instructString{3} = 'If smaller, press ''K'' (middle finger)';
    instructString{4} = 'Please respond as quickly as possible while still being accurate';
    instructString{5} = 'Press ''q'' to begin!';
    
    big=LEFT; small=RIGHT;
else
    instructString{1} = 'Decide whether each object is bigger or smaller than a shoe box in real life';
    instructString{2} = 'If smaller, press ''J'' (index finger)';
    instructString{3} = 'If bigger, press ''K'' (middle finger)';
    instructString{4} = 'Please respond as quickly as possible while still being accurate';
    instructString{5} = 'Press ''q'' to begin!';
    
    big=RIGHT; small=LEFT;
end


if practice==1
    instructString{end+1} = '(practice round)';
end

% instructions to get the experimenter at the end
if fMRI==0
    progressString{1} = 'Great job! Please find the';
    progressString{2} = 'experimenter in the other room.';
elseif fMRI==1
    progressString{1} = 'Great job!';
end

btw_repsString{1} = 'we''ll continue in a bit!';
btw_runs_time=7;%time to wait between runs
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

%% load trial sequences - based on the correct section
part='rem';
load([data_dir '/trial_sequences.mat'],'stim_seq_reminder','isi_reminder','stim_loc_seq_reminder','stim_in_pair_seq_reminder','reminder_onsets');
enc_idx=stim_seq_reminder;
stim_in_pair=stim_in_pair_seq_reminder;
loc_in_pair=stim_loc_seq_reminder;
trial_onsets=reminder_onsets;
isi_curr_part=isi_reminder;
%preparing a condition array
condition=rem(abs(enc_idx),10000);
condition=floor(condition/1000);

%preparing a stimuli size array
stim_size=floor(abs(enc_idx)/10000);

% if fMRI==0 || practice == 1
%     enc_trial_onsets=enc_trial_onsets-enc_trial_onsets(1,1,1)+run_delay; %#ok<NODEF>
% end

n_enc_runs=size(enc_idx,1);
%enc_idx is a matrixt with
%rows: number of sessions
%columns: order of presenting the stimuli - negative items should be responded
%with "smaller" and positivenumbers with "bigger"

if practice==1
    
    % override number of encoding trials and stim folder
    n_enc_runs=1;
    
    enc_idx=[1 -2 3 4 -5 6 -7 8];
    dir_list=dir([stim_dir '/practice']);
    
    % skip stuff that isn't an image
    dir_list = dir_list(3:end);               % skip . & ..
    if (strcmp(dir_list(1).name,'.DS_Store')==1) % also skip .DS_Store
        dir_list = dir_list(2:end);
    end
    
    prac_file_names={};
    % grab some objects
    for f=abs(enc_idx)
        prac_file_names{f}=[stim_dir '/practice/' dir_list(f).name];
    end
    
end

%preparing a correct_response array
correct_response=zeros(size(enc_idx));
correct_response(enc_idx>0)=big;
correct_response(enc_idx<0)=small;


%% load images
% create main window
mainWindow = Screen(windowIdx,'OpenWindow',backColor,[0 0 Screen_X Screen_Y]);
Screen(mainWindow,'TextFont',textFont);
Screen(mainWindow,'TextSize',textSize);

% stupid hack
if practice==1
    idx_to_find=unique(enc_idx);
    enc_texture=zeros(size(enc_idx));
    for i=idx_to_find %the number of the stimuli
        
        %enc_idx is the order of the stimuli
        % match encoding index to the list of all indices for this subject
        %all_idx is an array of numbers - possibly, this is the array that sets
        %which image is "image number x" it is different in each subject
        %takes "image number 1" finds it in loc - say, number 40.
        %temp_idx=all_idx==i;
        
        % read in images in the order of encoding presentation
        temp_image = imread(prac_file_names{abs(i)}); %now in temp image, the first 1, there is the 40th image in the folder
        
        % return index of images that can be called by 'DrawTexture'
        enc_texture(enc_idx==i) = Screen('MakeTexture',mainWindow,temp_image);
        %makes an array of images ready to be displayed
        %put in where there is the number "1" in enc_idx, the 40th image in the folder.
        %so, all_idx and he order in which items appear in enc_idx set the
        %randomization
    end
elseif practice==0
    idx_to_find=unique(enc_idx);
    enc_texture=zeros(size(enc_idx));
    for im=1:length(idx_to_find) %the number of the stimuli
        i=idx_to_find(im); %stupid hack because I don't understand why sometimes looping on arrays work and sometimes not.
        %enc_idx is the order of the stimuli
        % match encoding index to the list of all indices for this subject
        %all_idx is an array of numbers - possibly, this is the array that sets
        %which image is "image number x" it is different in each subject
        %takes "image number 1" finds it in loc - say, number 40.
        %temp_idx=all_idx==i;
        image_size=floor(abs(i)/10000);
        image_num=rem(abs(i),1000);
        % read in images in the order of encoding presentation
        temp_image = imread(file_names{image_num,image_size}); %now in temp image, the first 1, there is the 40th image in the folder
            
        % return index of images that can be called by 'DrawTexture'
        enc_texture(enc_idx==i) = Screen('MakeTexture',mainWindow,temp_image);
        %makes an array of images ready to be displayed
        %put in where there is the number "1" in enc_idx, the 40th image in the folder.
        %so, all_idx and he order in which items appear in enc_idx set the
        %randomization
    end
end

%% trial onsets and data structures

tpr=size(enc_texture,2); % trials per run

if debug==0 
    if practice==1 %this is practice - do not wait the initial 12 secs, just 2
        rand_onsets=ones(1,8)*0.5;
        rand_onsets=rand_onsets(randperm(length(rand_onsets)));
        trial_lengths=stim_dur + rand_onsets;
        trial_onsets=cumsum(trial_lengths,2) - (trial_lengths)+2;
        run_length=trial_onsets(end)+stim_dur+1; % end 1 sec after last trial
    end
elseif debug==1
    iti=.25 * ones(size(enc_texture));
    trial_lengths=stim_dur + iti;
    trial_onsets=cumsum(trial_lengths,2) - (stim_dur+iti);
    run_length=trial_onsets(end)+stim_dur+1; % end 1 sec after last trial
end


repeat =1;
while repeat==1
    
    stim_onset=zeros(n_enc_runs,tpr);
    response=zeros(n_enc_runs,tpr);
    response_times=zeros(n_enc_runs,tpr);
    response_accuracy=zeros(n_enc_runs,tpr);
    
    %% run loop
    for rep=startBlock:n_enc_runs
        if ~debug && ~practice %this is not practice, nor debugging
            run_length=trial_onsets(rep,end)+stim_dur+isi_curr_part(rep,end)+run_decay; % 10 secs added because scanner somehow too long
        end
        if practice == 0
            % set up output text file
            data_file=fullfile(data_dir,[sprintf('%s_',part) num2str(subj_num) subj_id '_run' num2str(rep) '.txt']);
            if exist(data_file,'file')
               data_file=[data_file(1:end-4) '1.txt'];
            end
            fid = fopen(data_file, 'w');
            fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', ...
                'subj', 'init', 'run', 'trial', 'onset', 'image', 'condition', 'pair', 'loc_in_pair', 'size', 'resp', 'acc', 'rt');
        end
        
        % print header to screen
        fprintf('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', ...
            'subj', 'init', 'run', 'trial', 'onset', 'image', 'condition', 'pair', 'loc_in_pair', 'size', 'resp', 'acc', 'rt');
        
        % set up diagnostic plot
        if fMRI==1 || debug==1
            rd = response_diagnostics(response_window, tpr);
        end
        
%         if practice==1
%             %% show PDF instructions
%             temp_image = imread('instructions/instructions_encoding.jpeg');
%             
%             % return index of images that can be called by 'DrawTexture'
%             instruct_texture = Screen('MakeTexture',mainWindow,temp_image);
%             
%             instruct_h=size(temp_image,1) ;
%             instruct_w=size(temp_image,2) ;
%             instructRect = [centerX-instruct_w/2,centerY-instruct_h/2,centerX+instruct_w/2,centerY+instruct_h/2];
%             
%             Screen('DrawTexture',mainWindow,instruct_texture,[0 0 instruct_w instruct_h],instructRect);
%             Screen('Flip',mainWindow);
%             
%             % wait for experimenter to advance with 'q' key
%             FlushEvents('keyDown');
%             while(1)
%                 temp = GetChar;
%                 if (temp == 'q')
%                     break;
%                 end
%             end
%         end
        
        %% show encoding instructions
        if rep==startBlock
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
            Priority(MaxPriority(windowIdx));
            
            if fMRI==1 && practice == 0
                
                % wait for scanner to spit out the first back tick
                FlushEvents('keyDown');
                while(1)
                    temp = GetChar;
                    if (temp == '`')
                        runTime=GetSecs;
                        break;
                    end
                end
                
            else
                runTime = GetSecs;
            end
        else
            runTime = GetSecs;
        end
        
        % fixation until first trial
        Screen(mainWindow,'FillOval',fixColor,fixDotRect);
        Screen('Flip',mainWindow);
        
        %% trial loop
        
        for trial=1:tpr
            
            % fixation
            Screen(mainWindow,'FillOval',fixColor,fixDotRect);
            Screen('Flip',mainWindow);
            
            % wait for next trial onset
            while GetSecs<runTime+trial_onsets(rep,trial); end
            
            % present images
            Screen('DrawTexture',mainWindow,enc_texture(rep,trial),imageRect,centerRect);
            stim_start=Screen('Flip',mainWindow);
            stim_onset(rep,trial)=stim_start-runTime;
            
            %  collect reponse
            removed=0;
            FlushEvents('keyDown');
            while GetSecs < stim_start + response_window
                WaitSecs(.0005);
                
                % take down image even if no response yet
                if ~removed && GetSecs >= stim_start + stim_dur
                    removed = 1;
                    Screen(mainWindow,'FillOval',fixColor,fixDotRect);
                    Screen('Flip', mainWindow);
                end
                
                % wait for response
                if (response(rep,trial) == 0)
                    
                    [keyIsDown, secs, keyCode] = KbCheck(device);
                    if (keyIsDown) && find(keyCode,1) ~= backtick
                        
                        % record response button and response time
                        response(rep,trial) = find(keyCode,1);
                        response_times(rep,trial) = secs - stim_start;
                        
                        % record accuracy
                        if find(keyCode,1) == correct_response(rep,trial)
                            response_accuracy(rep,trial) = 1;
                        else
                            response_accuracy(rep,trial) = 0;
                        end
                        
                    end
                end
            end
            
            % fixation until next trial
            Screen(mainWindow,'FillOval',fixColor,fixDotRect);
            Screen('Flip',mainWindow);
            
            % print output to screen and text file
            fprintf('%d\t%s\t%d\t%d\t%.3f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.3f\n', ...
                subj_num, subj_id, rep, trial, stim_onset(rep,trial), enc_idx(rep,trial), ...
                condition(rep,trial),stim_in_pair(rep,trial),loc_in_pair(rep,trial),stim_size(rep,trial),...
                response(rep,trial), response_accuracy(rep,trial), response_times(rep,trial)); %'condition', 'pair', 'loc_in_pair', 'size',
            
            if practice ==0
                fprintf(fid,'%d\t%s\t%d\t%d\t%.3f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.3f\n', ...
                    subj_num, subj_id, rep, trial, stim_onset(rep,trial), enc_idx(rep,trial), ...
                    condition(rep,trial),stim_in_pair(rep,trial),loc_in_pair(rep,trial),stim_size(rep,trial),...
                    response(rep,trial), response_accuracy(rep,trial), response_times(rep,trial));
            end
            
            if fMRI==1 || debug==1
                % plot performance and rt ongoing
                rd.log_trial(response_times(rep,trial), correct_response(rep,trial), response(rep,trial));
            end
        end
        
        % wait for scanner to finish
        while GetSecs < runTime+run_length; end
        
        
        if practice==0 
            % close this run's text file
            fclose(fid);
        end
        
        %have a short break if reps 1
        if rep == 1 && ~practice
            tempBounds = Screen('TextBounds',mainWindow,btw_repsString{1});
            Screen('drawtext',mainWindow,btw_repsString{1},round(centerX-tempBounds(3)/2),round(centerY-tempBounds(4)/5),textColor);
            clear tempBounds;
            Screen('Flip',mainWindow);
            while GetSecs < runTime+run_length+btw_runs_time; end
        end
        
        % repeat stimulus presentation during practice only if accuracy is not above 60%
        if practice==0 || (practice==1 && mean(response_accuracy(2:end))>.6)
            repeat=0;
        end
        
    end
end

%% clean up and go home
if  practice ==0
    mat_file=fullfile(data_dir,[sprintf('%s_',part) num2str(subj_num) subj_id '.mat']);
    if exist(mat_file,'file')
       mat_file=[mat_file(1:end-4) '1.mat'];
    end
    save(mat_file,'response','response_times','response_accuracy','stim_onset','description','condition','stim_in_pair','loc_in_pair','stim_size')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% now the error part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if practice==0
    
part='er_block';
load([data_dir '/trial_sequences.mat'],'stim_seq_vnv','isi_vnv_block','stim_loc_seq_vnv','stim_in_pair_seq_vnv','vnv_onsets');
enc_idx=stim_seq_vnv;
stim_in_pair=stim_in_pair_seq_vnv;
loc_in_pair=stim_loc_seq_vnv;
trial_onsets=vnv_onsets;
isi_curr_part=isi_vnv_block;
%preparing a condition array
condition=rem(abs(enc_idx),10000);
condition=floor(condition/1000);

%preparing a stimuli size array
stim_size=floor(abs(enc_idx)/10000);

% if fMRI==0 || practice == 1
%     enc_trial_onsets=enc_trial_onsets-enc_trial_onsets(1,1,1)+run_delay; %#ok<NODEF>
% end

n_enc_runs=size(enc_idx,1);
%enc_idx is a matrixt with
%rows: number of sessions
%columns: order of presenting the stimuli - negative items should be responded
%with "smaller" and positivenumbers with "bigger"

%preparing a correct_response array
correct_response=zeros(size(enc_idx));
correct_response(enc_idx>0)=big;
correct_response(enc_idx<0)=small;


%% load images
% create main window
mainWindow = Screen(windowIdx,'OpenWindow',backColor,[0 0 Screen_X Screen_Y]);
Screen(mainWindow,'TextFont',textFont);
Screen(mainWindow,'TextSize',textSize);

% stupid hack
idx_to_find=unique(enc_idx);
enc_texture=zeros(size(enc_idx));
for im=1:length(idx_to_find) %the number of the stimuli
    i=idx_to_find(im);
    %enc_idx is the order of the stimuli
    % match encoding index to the list of all indices for this subject
    %all_idx is an array of numbers - possibly, this is the array that sets
    %which image is "image number x" it is different in each subject
    %takes "image number 1" finds it in loc - say, number 40.
    %temp_idx=all_idx==i;
    image_size=floor(abs(i)/10000);
    image_num=rem(abs(i),1000);
    curr_cond=rem(abs(i),10000);
    curr_cond=floor(curr_cond/1000);
    % read in images in the order of encoding presentation
    if ismember(curr_cond,[5,6,7,8]) %take from the pairs
        temp_image = imread(file_names{image_num,image_size+numel(categ_names)});
    else
        temp_image = imread(file_names{image_num,image_size});
    end

    % return index of images that can be called by 'DrawTexture'
    enc_texture(enc_idx==i) = Screen('MakeTexture',mainWindow,temp_image);
    %makes an array of images ready to be displayed
    %put in where there is the number "1" in enc_idx, the 40th image in the folder.
    %so, all_idx and he order in which items appear in enc_idx set the
    %randomization
end

%% trial onsets and data structures

tpr=size(enc_texture,2); % trials per run

if debug==0 
    if practice==1 %this is practice - do not wait the initial 12 secs, just 2
        rand_onsets=[1 1 3 5 1 1 3 5];
        rand_onsets=rand_onsets(randperm(length(rand_onsets)));
        trial_lengths=stim_dur + rand_onsets;
        trial_onsets=cumsum(trial_lengths,2) - (trial_lengths)+2;
        run_length=trial_onsets(end)+stim_dur+1; % end 1 sec after last trial
    end
elseif debug==1
    iti=.25 * ones(size(enc_texture));
    trial_lengths=stim_dur + iti;
    trial_onsets=cumsum(trial_lengths,2) - (stim_dur+iti);
    run_length=trial_onsets(end)+stim_dur+1; % end 1 sec after last trial
end


repeat =1;
while repeat==1
    
    stim_onset=zeros(n_enc_runs,tpr);
    response=zeros(n_enc_runs,tpr);
    response_times=zeros(n_enc_runs,tpr);
    response_accuracy=zeros(n_enc_runs,tpr);
    
    %% run loop
    for rep=startBlock:n_enc_runs
        if ~debug && ~practice %this is not practice, nor debugging
            run_length=trial_onsets(rep,end)+stim_dur+isi_curr_part(rep,end)+run_decay; % 10 secs added because scanner somehow too long
        end
        
        % set up output text file
        data_file=fullfile(data_dir,[sprintf('%s_',part) num2str(subj_num) subj_id '_run' num2str(rep) '.txt']);
        if exist(data_file,'file')
            data_file=[data_file(1:end-4) '1.txt'];
        end
        fid = fopen(data_file, 'w');
        fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', ...
            'subj', 'init', 'run', 'trial', 'onset', 'image', 'condition', 'pair', 'loc_in_pair', 'size', 'resp', 'acc', 'rt');
        
        
        % print header to screen
        fprintf('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', ...
            'subj', 'init', 'run', 'trial', 'onset', 'image', 'condition', 'pair', 'loc_in_pair', 'size', 'resp', 'acc', 'rt');
        
        % set up diagnostic plot
        if fMRI==1 || debug==1
            rd = response_diagnostics(response_window, tpr);
        end
        
        
        %% show encoding instructions
        if rep==startBlock
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
            Priority(MaxPriority(windowIdx));
            
            if fMRI==1 && practice == 0
                
                % wait for scanner to spit out the first back tick
                FlushEvents('keyDown');
                while(1)
                    temp = GetChar;
                    if (temp == '`')
                        runTime=GetSecs;
                        break;
                    end
                end
                
            else
                runTime = GetSecs;
            end
        else
            runTime = GetSecs;
        end
        
        % fixation until first trial
        Screen(mainWindow,'FillOval',fixColor,fixDotRect);
        Screen('Flip',mainWindow);
        
        %% trial loop
        
        for trial=1:tpr
            
            % fixation
            Screen(mainWindow,'FillOval',fixColor,fixDotRect);
            Screen('Flip',mainWindow);
            
            % wait for next trial onset
            while GetSecs<runTime+trial_onsets(rep,trial); end
            
            % present images
            Screen('DrawTexture',mainWindow,enc_texture(rep,trial),imageRect,centerRect);
            stim_start=Screen('Flip',mainWindow);
            stim_onset(rep,trial)=stim_start-runTime;
            
            %  collect reponse
            removed=0;
            FlushEvents('keyDown');
            while GetSecs < stim_start + response_window
                WaitSecs(.0005);
                
                % take down image even if no response yet
                if ~removed && GetSecs >= stim_start + stim_dur
                    removed = 1;
                    Screen(mainWindow,'FillOval',fixColor,fixDotRect);
                    Screen('Flip', mainWindow);
                end
                
                % wait for response
                if (response(rep,trial) == 0)
                    
                    [keyIsDown, secs, keyCode] = KbCheck(device);
                    if (keyIsDown) && find(keyCode,1) ~= backtick
                        
                        % record response button and response time
                        response(rep,trial) = find(keyCode,1);
                        response_times(rep,trial) = secs - stim_start;
                        
                        % record accuracy
                        if find(keyCode,1) == correct_response(rep,trial)
                            response_accuracy(rep,trial) = 1;
                        else
                            response_accuracy(rep,trial) = 0;
                        end
                        
                    end
                end
            end
            
            % fixation until next trial
            Screen(mainWindow,'FillOval',fixColor,fixDotRect);
            Screen('Flip',mainWindow);
            
            % print output to screen and text file
            fprintf('%d\t%s\t%d\t%d\t%.3f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.3f\n', ...
                subj_num, subj_id, rep, trial, stim_onset(rep,trial), enc_idx(rep,trial), ...
                condition(rep,trial),stim_in_pair(rep,trial),loc_in_pair(rep,trial),stim_size(rep,trial),...
                response(rep,trial), response_accuracy(rep,trial), response_times(rep,trial)); %'condition', 'pair', 'loc_in_pair', 'size',
            
            if practice ==0
                fprintf(fid,'%d\t%s\t%d\t%d\t%.3f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.3f\n', ...
                    subj_num, subj_id, rep, trial, stim_onset(rep,trial), enc_idx(rep,trial), ...
                    condition(rep,trial),stim_in_pair(rep,trial),loc_in_pair(rep,trial),stim_size(rep,trial),...
                    response(rep,trial), response_accuracy(rep,trial), response_times(rep,trial));
            end
            
            if fMRI==1 || debug==1
                % plot performance and rt ongoing
                rd.log_trial(response_times(rep,trial), correct_response(rep,trial), response(rep,trial));
            end
        end
        
        % wait for scanner to finish
        while GetSecs < runTime+run_length; end
        
        if practice==0 
            % close this run's text file
            fclose(fid);
        end
        
        %have a short break if reps 1
        if rep < n_enc_runs && ~practice
            tempBounds = Screen('TextBounds',mainWindow,btw_repsString{1});
            Screen('drawtext',mainWindow,btw_repsString{1},round(centerX-tempBounds(3)/2),round(centerY-tempBounds(4)/5),textColor);
            clear tempBounds;
            Screen('Flip',mainWindow);
            while GetSecs < runTime+run_length+btw_runs_time; end
        end
        
        
        % repeat stimulus presentation during practice only if accuracy is not above 60%
        if practice==0 || (practice==1 && mean(response_accuracy(2:end))>.6)
            repeat=0;
        end
        
    end
end

end
%% clean up and go home
if  practice==0
    mat_file=fullfile(data_dir,[sprintf('%s_',part) num2str(subj_num) subj_id '.mat']);
    if exist(mat_file,'file')
       mat_file=[mat_file(1:end-4) '1.mat'];
    end
    save(mat_file,'response','response_times','response_accuracy','stim_onset','description','condition','stim_in_pair','loc_in_pair','stim_size')
end

if practice==0
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