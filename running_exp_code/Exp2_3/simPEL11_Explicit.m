function simPEL11_Explicit(subj_num, subj_id, varargin)
% function retrieval(subj_num, subj_initials, delay, fMRI, debug)
% Inputs:
%   subjectNumber: start first subject as 1, used to counterbalance button
%                       box response mappings (integer)
%   subj_initials:      subject initials, used to write output file names (string)
%   fMRI:               [0 / 1] looks for scanner pulse/trigger if fMRI,
%                       restricts KbCheck to button box
%   debug:              [0 / 1] if debug, not whole screen

%% boiler plate

fMRI=varargin{1};
practice=varargin{2};
debug=varargin{3};

if ~debug
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
    
    %that's for the regular ones:
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

% trial timing
if fMRI==1 && practice == 0
    run_delay = 12;
    run_decay = 12;
else
    run_delay = 2;
    run_decay = 2;
end

%subjects will have a maximum of 10 sec to make each decision.
if debug
    stim_dur = 0.10;           %sec
    response_window = 0.10;    %sec
    iti=0.1; %difference between trials
else
    stim_dur = 10;           %sec
    response_window = 10;    %sec
    iti=1; %difference between trials
end


% run_type={'AB','CD'};
% n_runs=2;

if length(varargin)==4
    start_run=varargin{4};
else
    start_run=1;
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

backtick=KbName('`~');
min_conf=KbName('1!');
max_conf=KbName('6^');

%confidence scale:
instructString{1} = 'Decide which of the bottom objects followed the top object during the bigger/smaller task';
instructString{2} = 'First choose one of the three bottom objects';
instructString{3} = 'Then rate your confidence 1-6';
instructString{4} = 'Please spread your ratings across all the 1-6 range';
instructString{5} = 'Press ''q'' to begin!';

% instructions to get the experimenter at the end
if fMRI==0
    progressString{1} = 'Great job! Please find the';
    progressString{2} = 'experimenter in the other room.';
elseif fMRI==1
    progressString{1} = 'Great job!';
end

%confidence scale:
% for i=1:6
%     confInst{i}=i;
% end
confInst{1} = 'very'; 
confInst{2} = 'unsure';
confInst{3} = 'sure';
confInst{4} = ' ';

%% set screen information
backColor = 127;
textColor = 0;
textColorFade = 64; %64 darker than screen, if want brighter: 224
textFont = 'Arial';
textSize = 18;
textSpacing = 25;
shift_up=30;%place figures slightly above the top/bottom idx to have space for the confedence scale
btm_pic_gap=50;%set the gap between the three bottom pics
conf_scale_bound=80;%set the bounds of the confidence scale, in pixels
fixColor = 0;
imageSize=350; % assumed square
%textSize=350; % assumed square
textImage_dist=10; %pixels
fixationSize = 4; % pixels
screens = Screen('Screens');

if fMRI
    % present stimuli in second screen
    windowIdx=max(screens);
    [Screen_X, Screen_Y]=Screen('WindowSize',windowIdx);
elseif debug==1
    % present stimuli on second screen, but smaller
    windowIdx=max(screens);
    Screen_X = 640;
    Screen_Y = 480;
elseif debug==2
    % present stimuli on second screen, but smaller
    windowIdx=max(screens);
    Screen_X = 1280;
    Screen_Y = 1024;
else
    % behavioral: run in testing room
    windowIdx=min(screens);
    [Screen_X, Screen_Y]=Screen('WindowSize',windowIdx);
end
centerX = Screen_X/2;
centerY = Screen_Y/2-shift_up;

rightX = Screen_X/2+btm_pic_gap+imageSize;
leftX  = Screen_X/2-btm_pic_gap-imageSize;
topY = Screen_Y/4-shift_up;
bottomY = Screen_Y*3/4-(shift_up*2);

% placeholder for images
imageRect = [0,0,imageSize,imageSize]; %x1,y1,x2,y2

% position of images
topRect = [centerX-imageSize/2,topY-imageSize/2,centerX+imageSize/2,topY+imageSize/2];
leftRect = [leftX-imageSize/2,bottomY-imageSize/2,leftX+imageSize/2,bottomY+imageSize/2];
rightRect = [rightX-imageSize/2,bottomY-imageSize/2,rightX+imageSize/2,bottomY+imageSize/2];
middleRect = [centerX-imageSize/2,bottomY-imageSize/2,centerX+imageSize/2,bottomY+imageSize/2];

%leftText= [[leftX-textSize/2,bottomY+imageSize/2+textImage_dist,leftX+textSize/2,bottomY+imageSize/2+textImage_dist+textSize];
% position of fixation dot
fixDotRect = [centerX-fixationSize,centerY-fixationSize,centerX+fixationSize,centerY+fixationSize];

%% load images - PART 1

% load priming trial sequences
load([data_dir '/Explicit_trial_sequencesPart1.mat'],'ExpTest_all_trials_InLoc','loc_allCondsInOrder','CondInOrder');

% create main window
mainWindow = Screen('OpenWindow',windowIdx,backColor,[0 0 Screen_X Screen_Y]);
Screen(mainWindow,'TextFont',textFont);
Screen(mainWindow,'TextSize',textSize);

Exp_texture=zeros(size(ExpTest_all_trials_InLoc));
Exp_idx=unique(ExpTest_all_trials_InLoc);
for i=1:length(Exp_idx) %for each row in ExpTest_all_trials_InLoc
    
    image_size=floor(abs(Exp_idx(i))/10000);
    image_num=rem(abs(Exp_idx(i)),1000);
    % read in images in the order of encoding presentation
    temp_image = imread(file_names{image_num,image_size}); %now in temp image, the first 1, there is the 40th image in the folder
    
    curr_cond=rem(abs(Exp_idx(i)),10000);
    curr_cond=floor(curr_cond/1000);
    % read in images in the order of encoding presentation
    if ismember(curr_cond,[5,6]) %take from the pairs
        temp_image = imread(file_names{image_num,image_size+numel(categ_names)});
    else
        temp_image = imread(file_names{image_num,image_size});
    end
    
    % return index of images that can be called by 'DrawTexture'
    Exp_texture(ExpTest_all_trials_InLoc==Exp_idx(i)) = Screen('MakeTexture',mainWindow,temp_image);
    
end

%% practice
if practice==1
    
    % override number of encoding trials and stim folder
    n_runs=1;
    prac_idx=[1 -2 3 -4 -5 6 7 8];
    dir_list=dir([stim_dir '/practice']);
    
    % skip stuff that isn't an image
    dir_list = dir_list(3:end);               % skip . & ..
    if (strcmp(dir_list(1).name,'.DS_Store')==1) % also skip .DS_Store
        dir_list = dir_list(2:end);
    end
    
    prac_file_names={};
    % grab some objects
    for f=abs(prac_idx)
        prac_file_names{f}=[stim_dir '/practice/' dir_list(f).name];
    end
    
    Exp_pract_inLoc=[-2 3 -4 1;6 7 8 -5];
    Exp_texture=zeros(size(Exp_pract_inLoc));
    Exp_idx=unique(Exp_pract_inLoc);
    for i=1:length(Exp_idx) %for each row in ExpTest_all_trials_InLoc
    curr_im=Exp_idx(i);    
    % read in images in the order of encoding presentation
    temp_image = imread(prac_file_names{abs(curr_im)}); %now in temp image, the first 1, there is the 40th image in the folder
    
    % return index of images that can be called by 'DrawTexture'
    Exp_texture(Exp_pract_inLoc==Exp_idx(i)) = Screen('MakeTexture',mainWindow,temp_image);
    
    end
    
end

%% trial onsets and data structures
%stupid hack: ALexa had 2 runs, I have one, just a way to get it without
%changing the script too much
start_run=1;
n_runs=1;

n_trials=size(Exp_texture,1);

stim_onset=zeros(n_trials,1);
conf_bt=zeros(n_trials,3); %the confidence response for each item
choice_bt=zeros(n_trials,3);%for each location
choice_rt=zeros(n_trials,3);%for each item - when did they choose it - from begning of the trial
conf_rt=zeros(n_trials,3);%for each item - duration of confidence rate, from time of choosing the item
curr_trial_dist=zeros(n_trials,3);%for each item, the distructors not based on their locations, but in the order of: 'C-target','R-tar/C-dist','C-dist',
curr_trial_loc=zeros(n_trials,3);%for each item, the location of the distructors in the order of: 'C-target','R-tar/C-dist','C-dist',
condition=zeros(n_trials,1);

for run=start_run:n_runs
    % set up output text file
    if practice==0
    data_file=fullfile(data_dir,['ExpTestPart1_' num2str(subj_num) subj_id '.txt']);
    if exist(data_file,'file')
       data_file=[data_file(1:end-4) '1.txt'];
    end
    fid = fopen(data_file, 'w');
    
    fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', ...
        'subj', 'init', 'trial', 'condition','onset','cue image','C-target','R-tar/C-dist','C-dist','loc: C-tar','loc: R-tar/C-dist','loc: C-dist',...
        'conf: target','conf: R-tar/C-dist','conf: C-dist','rt-conf: target','rt-conf: R-tar/C-dist','rt-conf: C-dist','bt-choice: target','bt-choice: R-tar/C-dist','bt-choice: C-dist','rt-choice: target','rt-choice: R-tar/C-dist','rt-choice: C-dist');
    
    % print header to screen
    fprintf('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', ...
        'subj', 'init', 'trial', 'condition','onset', 'cue image','C-target','R-tar/C-dist','C-dist','loc: C-tar','loc: R-tar/C-dist','loc: C-dist',...
        'conf: target','conf: R-tar/C-dist','conf: C-dist','rt-conf: target','rt-conf: R-tar/C-dist','rt-conf: C-dist','bt-choice: target','bt-choice: R-tar/C-dist','bt-choice: C-dist','rt-choice: target','rt-choice: R-tar/C-dist','rt-choice: C-dist');
    end
    % set up diagnostic plot
    if fMRI==1 || debug==1
        rd = response_diagnostics(response_window, n_trials);
    end
    

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
    
    % fixation until first trial
    fix_start=GetSecs;
    while(GetSecs< fix_start + iti)
        WaitSecs(.0005);
        Screen(mainWindow,'FillOval',fixColor,fixDotRect);
        Screen('Flip',mainWindow);
    end
    
    %% trial loop
    
    for trial=1:n_trials
        
        % fixation
        Screen(mainWindow,'FillOval',fixColor,fixDotRect);
        Screen('Flip',mainWindow);
        
        % wait for next trial onset
        %while GetSecs<runTime+trial_onsets(run,trial); end
        
        % present images
        Screen('DrawTexture',mainWindow,Exp_texture(trial,4),imageRect,topRect);
        Screen('DrawTexture',mainWindow,Exp_texture(trial,1),imageRect,leftRect);
        Screen('DrawTexture',mainWindow,Exp_texture(trial,2),imageRect,middleRect);
        Screen('DrawTexture',mainWindow,Exp_texture(trial,3),imageRect,rightRect);
        
        %current distractors (C-target, R-tar/C-dist, C-dist)
        curr_trial_dist(trial,:)=[ExpTest_all_trials_InLoc(trial,find(loc_allCondsInOrder(trial,:)==1)) ExpTest_all_trials_InLoc(trial,find(loc_allCondsInOrder(trial,:)==2)) ExpTest_all_trials_InLoc(trial,find(loc_allCondsInOrder(trial,:)==3))];
        %current locations (C-target, R-tar/C-dist, C-dist)
        curr_trial_loc(trial,:)=[find(loc_allCondsInOrder(trial,:)==1) find(loc_allCondsInOrder(trial,:)==2) find(loc_allCondsInOrder(trial,:)==3)];
        %current condition:
        condition(trial)=rem(abs(ExpTest_all_trials_InLoc(trial,4)),10000);
        condition=floor(abs(condition/1000));
        
        stim_start=Screen('Flip',mainWindow,0,1);
        stim_onset(trial)=stim_start-runTime;
        resp_start=stim_start;%counter for the time for each response
        
        %  collect reponse
        removed=0;
        curr_trial_responses=0;% a slot for each response
        curr_responded_stim=1;%first stimuli to get response, update when a response is made
        
        FlushEvents('keyDown');
        while (GetSecs < resp_start + response_window) && removed==0
            WaitSecs(.0005);
            
            % take down image even if no response yet
            if ~removed && GetSecs >= resp_start + response_window
                removed = 1;
                Screen(mainWindow,'FillOval',fixColor,fixDotRect);
                Screen('Flip', mainWindow);
            end
            
            
            % wait for response
            if (curr_trial_responses(curr_responded_stim) == 0) && (~removed) %no selection of a distractor was made
                
                [keyIsDown, secs, keyCode] = KbCheck(device);
                if (keyIsDown) && find(keyCode,1) ~= backtick
                    curr_key=find(keyCode,1);
                    
                    switch curr_key
                        case LEFT %subject chose the left one
                            %make sure that they haven't mistakenly pressed
                            %a chosen dist:
                            if any(curr_trial_responses==1)
                                FlushEvents('keyDown');
                                clear keyIsDown secs keyCode;
                            else
                                curr_trial_responses(curr_responded_stim)=1;
                                curr_dist=loc_allCondsInOrder(trial,1); %can be 1/2/3
                                center_conf_text=leftX; %set for the location of the confidence instructions
                                % record response button and response time
                                choice_bt(trial,curr_dist) = find(keyCode,1);
                                choice_rt(trial,curr_dist) = secs - stim_start; %from begining of the trial
                            end
                        case MIDDLE
                            if any(curr_trial_responses==2)
                                FlushEvents('keyDown');
                                clear keyIsDown secs keyCode;
                            else
                                curr_trial_responses(curr_responded_stim)=2;
                                curr_dist=loc_allCondsInOrder(trial,2); %can be 1/2/3
                                center_conf_text=centerX; %set for the location of the confidence instructions
                                % record response button and response time
                                choice_bt(trial,curr_dist) = find(keyCode,1);
                                choice_rt(trial,curr_dist) = secs - stim_start; %from begining of the trial
                            end
                        case RIGHT
                            if any(curr_trial_responses==3)
                                FlushEvents('keyDown');
                                clear keyIsDown secs keyCode;
                            else
                                curr_trial_responses(curr_responded_stim)=3;
                                curr_dist=loc_allCondsInOrder(trial,3); %can be 1/2/3
                                center_conf_text=rightX; %set for the location of the confidence instructions
                                % record response button and response time
                                choice_bt(trial,curr_dist) = find(keyCode,1);
                                choice_rt(trial,curr_dist) = secs - stim_start; %from begining of the trial
                            end
                    end
                    
                    
                    if (curr_trial_responses(curr_responded_stim) ~= 0) %subject chose a distructor
                        
                        %update that a response was given
                        curr_responded_stim=curr_responded_stim+1;
                        
                        %set the confidence text
                        for i=1:6
                            Screen('drawtext',mainWindow,num2str(i),round(center_conf_text-(imageSize-conf_scale_bound)/2+((i-1)*((imageSize-conf_scale_bound)/5))-5),round(bottomY+imageSize/2+textImage_dist+textSpacing),textColor);
                        end
                        %
                        
                        tempBounds = Screen('TextBounds',mainWindow,confInst{4}); %followed?
                        Screen('drawtext',mainWindow,confInst{4},round(center_conf_text-tempBounds(3)/2),round(bottomY+imageSize/2+textImage_dist-tempBounds(4)/5),textColor);
                        clear tempBounds;
                        
                        tempBounds = Screen('TextBounds',mainWindow,confInst{1});
                        Screen('drawtext',mainWindow,confInst{1},round(center_conf_text-(imageSize-conf_scale_bound)/2-tempBounds(3)/2),round(bottomY+imageSize/2+textImage_dist-tempBounds(4)/5+(2*textSpacing)),textColor);
                        clear tempBounds;
                        
                        tempBounds = Screen('TextBounds',mainWindow,confInst{2});
                        Screen('drawtext',mainWindow,confInst{2},round(center_conf_text-(imageSize-conf_scale_bound)/2-tempBounds(3)/2),round(bottomY+imageSize/2+textImage_dist-tempBounds(4)/5+(3*textSpacing-3)),textColor);
                        clear tempBounds;
                        
                        tempBounds = Screen('TextBounds',mainWindow,confInst{1});
                        Screen('drawtext',mainWindow,confInst{1},round(center_conf_text+(imageSize-conf_scale_bound)/2-tempBounds(3)/2),round(bottomY+imageSize/2+textImage_dist-tempBounds(4)/5+(2*textSpacing)),textColor);
                        clear tempBounds;
                        
                        tempBounds = Screen('TextBounds',mainWindow,confInst{3});
                        Screen('drawtext',mainWindow,confInst{3},round(center_conf_text+(imageSize-conf_scale_bound)/2-tempBounds(3)/2),round(bottomY+imageSize/2+textImage_dist-tempBounds(4)/5+(3*textSpacing-3)),textColor);
                        clear tempBounds;
                        
                        resp_start=Screen('Flip', mainWindow,0,1);
                        
                        %wait before continue - stuipd hack so that ptb
                        %won't count the response of choosing the stimuli
                        %as the response to the confidence
                        while GetSecs<resp_start+0.2; end 
                        
                        
                        gave_conf_resp=0;
                        FlushEvents('keyDown');
                        clear keyIsDown secs keyCode;
                        while (GetSecs < (resp_start + response_window)) && gave_conf_resp==0
                            WaitSecs(.0005);
                            %entered_conf=GetSecs
                            % fade confidence scale even if no response yet
                            if GetSecs >= (resp_start + stim_dur)
                                %entered_too_long=GetSecs
                                %set a grey rectengle to hide the conf scale before
                                %repainting it in grey:
                                %grap tempBounds from the instructions to know how
                                %big it should be:
                                tempBounds = Screen('TextBounds',mainWindow,confInst{3});
                                Screen('FillRect',mainWindow,backColor,[round(center_conf_text-(imageSize-conf_scale_bound)/2-tempBounds(3)/2)-2,round(bottomY+imageSize/2),round(center_conf_text+(imageSize-conf_scale_bound)/2+tempBounds(3)/2+2),round(bottomY+imageSize/2+textImage_dist+tempBounds(4)/2+(3*textSpacing+3)+2)]);
                                clear tempBounds;
                                
                                %grey the current confidence text
                                %set the confidence text
                                for i=1:6
                                    Screen('drawtext',mainWindow,num2str(i),round(center_conf_text-(imageSize-conf_scale_bound)/2+((i-1)*((imageSize-conf_scale_bound)/5))-5),round(bottomY+imageSize/2+textImage_dist+textSpacing),textColorFade);
                                end
                                %
                                tempBounds = Screen('TextBounds',mainWindow,confInst{4}); %followed?
                                Screen('drawtext',mainWindow,confInst{4},round(center_conf_text-tempBounds(3)/2),round(bottomY+imageSize/2+textImage_dist-tempBounds(4)/5),textColorFade);
                                clear tempBounds;
                                
                                tempBounds = Screen('TextBounds',mainWindow,confInst{1});
                                Screen('drawtext',mainWindow,confInst{1},round(center_conf_text-(imageSize-conf_scale_bound)/2-tempBounds(3)/2),round(bottomY+imageSize/2+textImage_dist-tempBounds(4)/5+(2*textSpacing)),textColorFade);
                                clear tempBounds;
                                
                                tempBounds = Screen('TextBounds',mainWindow,confInst{2});
                                Screen('drawtext',mainWindow,confInst{2},round(center_conf_text-(imageSize-conf_scale_bound)/2-tempBounds(3)/2),round(bottomY+imageSize/2+textImage_dist-tempBounds(4)/5+(3*textSpacing-3)),textColorFade);
                                clear tempBounds;
                                
                                tempBounds = Screen('TextBounds',mainWindow,confInst{1});
                                Screen('drawtext',mainWindow,confInst{1},round(center_conf_text+(imageSize-conf_scale_bound)/2-tempBounds(3)/2),round(bottomY+imageSize/2+textImage_dist-tempBounds(4)/5+(2*textSpacing)),textColorFade);
                                clear tempBounds;
                                
                                tempBounds = Screen('TextBounds',mainWindow,confInst{3});
                                Screen('drawtext',mainWindow,confInst{3},round(center_conf_text+(imageSize-conf_scale_bound)/2-tempBounds(3)/2),round(bottomY+imageSize/2+textImage_dist-tempBounds(4)/5+(3*textSpacing-3)),textColorFade);
                                clear tempBounds;
                                
                                %turn it al grey and update response time
                                resp_start=Screen('Flip', mainWindow,0,1);
                                %wait before continue - stuipd hack so that ptb
                                %won't count the response of choosing the stimuli
                                %as the response to the confidence
                                while GetSecs<resp_start+0.1; end
                                
                                
                                %get out of the confidence loop
                                gave_conf_resp=1;
                                
                            end
                            
                            % wait for response
                            if (conf_bt(trial,curr_dist) == 0) && gave_conf_resp==0 %no confidence judgement was made for the current distructor
                                %disp('waiting conf resp')
                                [keyIsDown, secs, keyCode] = KbCheck(device);
                                if (keyIsDown) && (find(keyCode,1)-min_conf+1>=1 && find(keyCode,1)-max_conf+1<=6) %subject gave a confidence response
                                
                                    %disp('entered if for got conf resp')
                                    conf_bt(trial,curr_dist) = find(keyCode,1)-min_conf+1;%the codes of the numbers start from 89=1, 90-2.., tha will write the specific confidence rating in the logfile
                                    conf_rt(trial,curr_dist) = secs - resp_start; %from choosing the current dist;
                                    
                                    %set a grey rectengle to hide the conf scale before
                                    %repainting it in grey:
                                    %grap tempBounds from the instructions to know how
                                    %big it should be:
                                    tempBounds = Screen('TextBounds',mainWindow,confInst{3});
                                    Screen('FillRect',mainWindow,backColor,[round(center_conf_text-(imageSize-conf_scale_bound)/2-tempBounds(3)/2)-2,round(bottomY+imageSize/2),round(center_conf_text+(imageSize-conf_scale_bound)/2+tempBounds(3)/2+2),round(bottomY+imageSize/2+textImage_dist+tempBounds(4)/2+(3*textSpacing+3)+2)]);
                                    clear tempBounds;
                                    
                                    %grey the current confidence text
                                    %set the confidence text
                                    for i=1:6
                                        Screen('drawtext',mainWindow,num2str(i),round(center_conf_text-(imageSize-conf_scale_bound)/2+((i-1)*((imageSize-conf_scale_bound)/5))-5),round(bottomY+imageSize/2+textImage_dist+textSpacing),textColorFade);
                                    end
                                    %
                                    tempBounds = Screen('TextBounds',mainWindow,confInst{4}); %followed?
                                    Screen('drawtext',mainWindow,confInst{4},round(center_conf_text-tempBounds(3)/2),round(bottomY+imageSize/2+textImage_dist-tempBounds(4)/5),textColorFade);
                                    clear tempBounds;
                                    
                                    tempBounds = Screen('TextBounds',mainWindow,confInst{1});
                                    Screen('drawtext',mainWindow,confInst{1},round(center_conf_text-(imageSize-conf_scale_bound)/2-tempBounds(3)/2),round(bottomY+imageSize/2+textImage_dist-tempBounds(4)/5+(2*textSpacing)),textColorFade);
                                    clear tempBounds;
                                    
                                    tempBounds = Screen('TextBounds',mainWindow,confInst{2});
                                    Screen('drawtext',mainWindow,confInst{2},round(center_conf_text-(imageSize-conf_scale_bound)/2-tempBounds(3)/2),round(bottomY+imageSize/2+textImage_dist-tempBounds(4)/5+(3*textSpacing-3)),textColorFade);
                                    clear tempBounds;
                                    
                                    tempBounds = Screen('TextBounds',mainWindow,confInst{1});
                                    Screen('drawtext',mainWindow,confInst{1},round(center_conf_text+(imageSize-conf_scale_bound)/2-tempBounds(3)/2),round(bottomY+imageSize/2+textImage_dist-tempBounds(4)/5+(2*textSpacing)),textColorFade);
                                    clear tempBounds;
                                    
                                    tempBounds = Screen('TextBounds',mainWindow,confInst{3});
                                    Screen('drawtext',mainWindow,confInst{3},round(center_conf_text+(imageSize-conf_scale_bound)/2-tempBounds(3)/2),round(bottomY+imageSize/2+textImage_dist-tempBounds(4)/5+(3*textSpacing-3)),textColorFade);
                                    clear tempBounds;
                                    
                                    %turn it al grey and update response time
                                    resp_start=Screen('Flip', mainWindow,0,1);
                                    
                                    %wait before continue - stuipd hack so that ptb
                                    %won't count the response of choosing the stimuli
                                    %as the response to the confidence
                                    while GetSecs<resp_start+0.2; end
                                  
                                    %get out of the confidence loop
                                    gave_conf_resp=1;
                                    
                                end
                                
                            end %ends if a conf response was given
                            
                        end %ends the while loop for duration of the trial
                    end %ends the conditional of wether the subject gave a response for the distractor
                    
                end
            end %ends if a distractor was chosen
            
            if (curr_trial_responses~=0) %subject gave 3 responses - move to the next trial
                removed=1;
                %clear the screen
                Screen(mainWindow,'FillOval',fixColor,fixDotRect);
                Screen('Flip',mainWindow);
            end
        end %ends the while loop for duration of the trial
        
        % fixation until next trial
        fix_start=GetSecs;
        while(GetSecs< fix_start + iti)
            WaitSecs(.0005);
            Screen(mainWindow,'FillOval',fixColor,fixDotRect);
            Screen('Flip',mainWindow);
        end
        
        % print output to screen and text file
        if practice==0
        fprintf('%d\t%s\t%d\t%d\t%.3f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.3f\t%.3f\t%.3f\t%d\t%d\t%d\t%.3f\t%.3f\t%.3f\n', ...
            subj_num, subj_id, trial, condition(trial), stim_onset(trial), ExpTest_all_trials_InLoc(trial,4), curr_trial_dist(trial,:),curr_trial_loc(trial,:),...
            conf_bt(trial,:),conf_rt(trial,:),choice_bt(trial,:), choice_rt(trial,:));
        
        fprintf(fid,'%d\t%s\t%d\t%d\t%.3f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.3f\t%.3f\t%.3f\t%d\t%d\t%d\t%.3f\t%.3f\t%.3f\n', ...
            subj_num, subj_id, trial, condition(trial), stim_onset(trial), ExpTest_all_trials_InLoc(trial,4), curr_trial_dist(trial,:),curr_trial_loc(trial,:),...
            conf_bt(trial,:),conf_rt(trial,:),choice_bt(trial,:), choice_rt(trial,:));
        end
        % plot performance and rt ongoing
%         if fMRI==1 || debug==1
%             rd.log_trial(response_times(run,trial), correct_response(run,trial), response(run,trial));
%         end
    end
    
    % wait for scanner to finish
    %while GetSecs < run_length; end
    
    % close this run's text file
    if practice==0
        fclose(fid);
    end
end % run loop

%% clean up and go home
if practice==0
    mat_file=fullfile(data_dir,'explicitPart1.mat');
    if exist(mat_file,'file')
       mat_file=[mat_file(1:end-4) '1.mat'];
    end
    save(mat_file,'condition','stim_onset','ExpTest_all_trials_InLoc','curr_trial_dist','curr_trial_loc','conf_bt','conf_rt','choice_bt','choice_rt','description')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load images - PART 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load priming trial sequences
if practice==0
load([data_dir '/Explicit_trial_sequencesPart2.mat'],'ExpTest_all_trials_InLoc','loc_allCondsInOrder','CondInOrder');

Exp_texture=zeros(size(ExpTest_all_trials_InLoc));
Exp_idx=unique(ExpTest_all_trials_InLoc);
for i=1:length(Exp_idx) %for each row in ExpTest_all_trials_InLoc
    
    image_size=floor(abs(Exp_idx(i))/10000);
    image_num=rem(abs(Exp_idx(i)),1000);
    % read in images in the order of encoding presentation
    temp_image = imread(file_names{image_num,image_size}); %now in temp image, the first 1, there is the 40th image in the folder
    
    % return index of images that can be called by 'DrawTexture'
    Exp_texture(ExpTest_all_trials_InLoc==Exp_idx(i)) = Screen('MakeTexture',mainWindow,temp_image);
    
end

%% trial onsets and data structures
%stupid hack: ALexa had 2 runs, I have one, just a way to get it without
%changing the script too much
start_run=1;
n_runs=1;

n_trials=size(Exp_texture,1);

stim_onset=zeros(n_trials,1);
conf_bt=zeros(n_trials,3); %the confidence response for each item
choice_bt=zeros(n_trials,3);%for each location
choice_rt=zeros(n_trials,3);%for each item - when did they choose it - from begning of the trial
conf_rt=zeros(n_trials,3);%for each item - duration of confidence rate, from time of choosing the item
curr_trial_dist=zeros(n_trials,3);%for each item, the distructors not based on their locations, but in the order of: 'C-target','R-tar/C-dist','C-dist',
curr_trial_loc=zeros(n_trials,3);%for each item, the location of the distructors in the order of: 'C-target','R-tar/C-dist','C-dist',
condition=zeros(n_trials,1);

for run=start_run:n_runs
    % set up output text file
    if practice==0
    data_file=fullfile(data_dir,['ExpTestPart2_' num2str(subj_num) subj_id '.txt']);
    if exist(data_file,'file')
       data_file=[data_file(1:end-4) '1.txt'];
    end
    fid = fopen(data_file, 'w');
    
    fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', ...
        'subj', 'init', 'trial', 'condition','onset','cue image','C-target','R-tar/C-dist','C-dist','loc: C-tar','loc: R-tar/C-dist','loc: C-dist',...
        'conf: target','conf: R-tar/C-dist','conf: C-dist','rt-conf: target','rt-conf: R-tar/C-dist','rt-conf: C-dist','bt-choice: target','bt-choice: R-tar/C-dist','bt-choice: C-dist','rt-choice: target','rt-choice: R-tar/C-dist','rt-choice: C-dist');
    
    % print header to screen
    fprintf('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', ...
        'subj', 'init', 'trial', 'condition','onset', 'cue image','C-target','R-tar/C-dist','C-dist','loc: C-tar','loc: R-tar/C-dist','loc: C-dist',...
        'conf: target','conf: R-tar/C-dist','conf: C-dist','rt-conf: target','rt-conf: R-tar/C-dist','rt-conf: C-dist','bt-choice: target','bt-choice: R-tar/C-dist','bt-choice: C-dist','rt-choice: target','rt-choice: R-tar/C-dist','rt-choice: C-dist');
    end
    % set up diagnostic plot
    if fMRI==1 || debug==1
        rd = response_diagnostics(response_window, n_trials);
    end
    

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
    
    % fixation until first trial
    fix_start=GetSecs;
    while(GetSecs< fix_start + iti)
        WaitSecs(.0005);
        Screen(mainWindow,'FillOval',fixColor,fixDotRect);
        Screen('Flip',mainWindow);
    end
    
    %% trial loop
    
    for trial=1:n_trials
        
        % fixation
        Screen(mainWindow,'FillOval',fixColor,fixDotRect);
        Screen('Flip',mainWindow);
        
        % wait for next trial onset
        %while GetSecs<runTime+trial_onsets(run,trial); end
        
        % present images
        Screen('DrawTexture',mainWindow,Exp_texture(trial,4),imageRect,topRect);
        Screen('DrawTexture',mainWindow,Exp_texture(trial,1),imageRect,leftRect);
        Screen('DrawTexture',mainWindow,Exp_texture(trial,2),imageRect,middleRect);
        Screen('DrawTexture',mainWindow,Exp_texture(trial,3),imageRect,rightRect);
        
        %current distractors (C-target, R-tar/C-dist, C-dist)
        curr_trial_dist(trial,:)=[ExpTest_all_trials_InLoc(trial,find(loc_allCondsInOrder(trial,:)==1)) ExpTest_all_trials_InLoc(trial,find(loc_allCondsInOrder(trial,:)==2)) ExpTest_all_trials_InLoc(trial,find(loc_allCondsInOrder(trial,:)==3))];
        %current locations (C-target, R-tar/C-dist, C-dist)
        curr_trial_loc(trial,:)=[find(loc_allCondsInOrder(trial,:)==1) find(loc_allCondsInOrder(trial,:)==2) find(loc_allCondsInOrder(trial,:)==3)];
        %current condition:
        condition(trial)=rem(abs(ExpTest_all_trials_InLoc(trial,4)),10000);
        condition=floor(abs(condition/1000));
        
        stim_start=Screen('Flip',mainWindow,0,1);
        stim_onset(trial)=stim_start-runTime;
        resp_start=stim_start;%counter for the time for each response
        
        %  collect reponse
        removed=0;
        curr_trial_responses=0;% a slot for each response
        curr_responded_stim=1;%first stimuli to get response, update when a response is made
        
        FlushEvents('keyDown');
        while (GetSecs < resp_start + response_window) && removed==0
            WaitSecs(.0005);
            
            % take down image even if no response yet
            if ~removed && GetSecs >= resp_start + response_window
                removed = 1;
                Screen(mainWindow,'FillOval',fixColor,fixDotRect);
                Screen('Flip', mainWindow);
            end
            
            
            % wait for response
            if (curr_trial_responses(curr_responded_stim) == 0) && (~removed) %no selection of a distractor was made
                
                [keyIsDown, secs, keyCode] = KbCheck(device);
                if (keyIsDown) && find(keyCode,1) ~= backtick
                    curr_key=find(keyCode,1);
                    
                    switch curr_key
                        case LEFT %subject chose the left one
                            %make sure that they haven't mistakenly pressed
                            %a chosen dist:
                            if any(curr_trial_responses==1)
                                FlushEvents('keyDown');
                                clear keyIsDown secs keyCode;
                            else
                                curr_trial_responses(curr_responded_stim)=1;
                                curr_dist=loc_allCondsInOrder(trial,1); %can be 1/2/3
                                center_conf_text=leftX; %set for the location of the confidence instructions
                                % record response button and response time
                                choice_bt(trial,curr_dist) = find(keyCode,1);
                                choice_rt(trial,curr_dist) = secs - stim_start; %from begining of the trial
                            end
                        case MIDDLE
                            if any(curr_trial_responses==2)
                                FlushEvents('keyDown');
                                clear keyIsDown secs keyCode;
                            else
                                curr_trial_responses(curr_responded_stim)=2;
                                curr_dist=loc_allCondsInOrder(trial,2); %can be 1/2/3
                                center_conf_text=centerX; %set for the location of the confidence instructions
                                % record response button and response time
                                choice_bt(trial,curr_dist) = find(keyCode,1);
                                choice_rt(trial,curr_dist) = secs - stim_start; %from begining of the trial
                            end
                        case RIGHT
                            if any(curr_trial_responses==3)
                                FlushEvents('keyDown');
                                clear keyIsDown secs keyCode;
                            else
                                curr_trial_responses(curr_responded_stim)=3;
                                curr_dist=loc_allCondsInOrder(trial,3); %can be 1/2/3
                                center_conf_text=rightX; %set for the location of the confidence instructions
                                % record response button and response time
                                choice_bt(trial,curr_dist) = find(keyCode,1);
                                choice_rt(trial,curr_dist) = secs - stim_start; %from begining of the trial
                            end
                    end
                    
                    
                    if (curr_trial_responses(curr_responded_stim) ~= 0) %subject chose a distructor
                        
                        %update that a response was given
                        curr_responded_stim=curr_responded_stim+1;
                        
                        %set the confidence text
                        for i=1:6
                            Screen('drawtext',mainWindow,num2str(i),round(center_conf_text-(imageSize-conf_scale_bound)/2+((i-1)*((imageSize-conf_scale_bound)/5))-5),round(bottomY+imageSize/2+textImage_dist+textSpacing),textColor);
                        end
                        %
                        
                        tempBounds = Screen('TextBounds',mainWindow,confInst{4}); %followed?
                        Screen('drawtext',mainWindow,confInst{4},round(center_conf_text-tempBounds(3)/2),round(bottomY+imageSize/2+textImage_dist-tempBounds(4)/5),textColor);
                        clear tempBounds;
                        
                        tempBounds = Screen('TextBounds',mainWindow,confInst{1});
                        Screen('drawtext',mainWindow,confInst{1},round(center_conf_text-(imageSize-conf_scale_bound)/2-tempBounds(3)/2),round(bottomY+imageSize/2+textImage_dist-tempBounds(4)/5+(2*textSpacing)),textColor);
                        clear tempBounds;
                        
                        tempBounds = Screen('TextBounds',mainWindow,confInst{2});
                        Screen('drawtext',mainWindow,confInst{2},round(center_conf_text-(imageSize-conf_scale_bound)/2-tempBounds(3)/2),round(bottomY+imageSize/2+textImage_dist-tempBounds(4)/5+(3*textSpacing-3)),textColor);
                        clear tempBounds;
                        
                        tempBounds = Screen('TextBounds',mainWindow,confInst{1});
                        Screen('drawtext',mainWindow,confInst{1},round(center_conf_text+(imageSize-conf_scale_bound)/2-tempBounds(3)/2),round(bottomY+imageSize/2+textImage_dist-tempBounds(4)/5+(2*textSpacing)),textColor);
                        clear tempBounds;
                        
                        tempBounds = Screen('TextBounds',mainWindow,confInst{3});
                        Screen('drawtext',mainWindow,confInst{3},round(center_conf_text+(imageSize-conf_scale_bound)/2-tempBounds(3)/2),round(bottomY+imageSize/2+textImage_dist-tempBounds(4)/5+(3*textSpacing-3)),textColor);
                        clear tempBounds;
                        
                        resp_start=Screen('Flip', mainWindow,0,1);
                        
                        %wait before continue - stuipd hack so that ptb
                        %won't count the response of choosing the stimuli
                        %as the response to the confidence
                        while GetSecs<resp_start+0.2; end 
                        
                        
                        gave_conf_resp=0;
                        FlushEvents('keyDown');
                        clear keyIsDown secs keyCode;
                        while (GetSecs < (resp_start + response_window)) && gave_conf_resp==0
                            WaitSecs(.0005);
                            %entered_conf=GetSecs
                            % fade confidence scale even if no response yet
                            if GetSecs >= (resp_start + stim_dur)
                                %entered_too_long=GetSecs
                                %set a grey rectengle to hide the conf scale before
                                %repainting it in grey:
                                %grap tempBounds from the instructions to know how
                                %big it should be:
                                tempBounds = Screen('TextBounds',mainWindow,confInst{3});
                                Screen('FillRect',mainWindow,backColor,[round(center_conf_text-(imageSize-conf_scale_bound)/2-tempBounds(3)/2)-2,round(bottomY+imageSize/2),round(center_conf_text+(imageSize-conf_scale_bound)/2+tempBounds(3)/2+2),round(bottomY+imageSize/2+textImage_dist+tempBounds(4)/2+(3*textSpacing+3)+2)]);
                                clear tempBounds;
                                
                                %grey the current confidence text
                                %set the confidence text
                                for i=1:6
                                    Screen('drawtext',mainWindow,num2str(i),round(center_conf_text-(imageSize-conf_scale_bound)/2+((i-1)*((imageSize-conf_scale_bound)/5))-5),round(bottomY+imageSize/2+textImage_dist+textSpacing),textColorFade);
                                end
                                %
                                tempBounds = Screen('TextBounds',mainWindow,confInst{4}); %followed?
                                Screen('drawtext',mainWindow,confInst{4},round(center_conf_text-tempBounds(3)/2),round(bottomY+imageSize/2+textImage_dist-tempBounds(4)/5),textColorFade);
                                clear tempBounds;
                                
                                tempBounds = Screen('TextBounds',mainWindow,confInst{1});
                                Screen('drawtext',mainWindow,confInst{1},round(center_conf_text-(imageSize-conf_scale_bound)/2-tempBounds(3)/2),round(bottomY+imageSize/2+textImage_dist-tempBounds(4)/5+(2*textSpacing)),textColorFade);
                                clear tempBounds;
                                
                                tempBounds = Screen('TextBounds',mainWindow,confInst{2});
                                Screen('drawtext',mainWindow,confInst{2},round(center_conf_text-(imageSize-conf_scale_bound)/2-tempBounds(3)/2),round(bottomY+imageSize/2+textImage_dist-tempBounds(4)/5+(3*textSpacing-3)),textColorFade);
                                clear tempBounds;
                                
                                tempBounds = Screen('TextBounds',mainWindow,confInst{1});
                                Screen('drawtext',mainWindow,confInst{1},round(center_conf_text+(imageSize-conf_scale_bound)/2-tempBounds(3)/2),round(bottomY+imageSize/2+textImage_dist-tempBounds(4)/5+(2*textSpacing)),textColorFade);
                                clear tempBounds;
                                
                                tempBounds = Screen('TextBounds',mainWindow,confInst{3});
                                Screen('drawtext',mainWindow,confInst{3},round(center_conf_text+(imageSize-conf_scale_bound)/2-tempBounds(3)/2),round(bottomY+imageSize/2+textImage_dist-tempBounds(4)/5+(3*textSpacing-3)),textColorFade);
                                clear tempBounds;
                                
                                %turn it al grey and update response time
                                resp_start=Screen('Flip', mainWindow,0,1);
                                %wait before continue - stuipd hack so that ptb
                                %won't count the response of choosing the stimuli
                                %as the response to the confidence
                                while GetSecs<resp_start+0.1; end
                                
                                
                                %get out of the confidence loop
                                gave_conf_resp=1;
                                
                            end
                            
                            % wait for response
                            if (conf_bt(trial,curr_dist) == 0) && gave_conf_resp==0 %no confidence judgement was made for the current distructor
                                %disp('waiting conf resp')
                                [keyIsDown, secs, keyCode] = KbCheck(device);
                                if (keyIsDown) && (find(keyCode,1)-min_conf+1>=1 && find(keyCode,1)-max_conf+1<=6) %subject gave a confidence response
                                
                                    %disp('entered if for got conf resp')
                                    conf_bt(trial,curr_dist) = find(keyCode,1)-min_conf+1;%the codes of the numbers start from 89=1, 90-2.., tha will write the specific confidence rating in the logfile
                                    conf_rt(trial,curr_dist) = secs - resp_start; %from choosing the current dist;
                                    
                                    %set a grey rectengle to hide the conf scale before
                                    %repainting it in grey:
                                    %grap tempBounds from the instructions to know how
                                    %big it should be:
                                    tempBounds = Screen('TextBounds',mainWindow,confInst{3});
                                    Screen('FillRect',mainWindow,backColor,[round(center_conf_text-(imageSize-conf_scale_bound)/2-tempBounds(3)/2)-2,round(bottomY+imageSize/2),round(center_conf_text+(imageSize-conf_scale_bound)/2+tempBounds(3)/2+2),round(bottomY+imageSize/2+textImage_dist+tempBounds(4)/2+(3*textSpacing+3)+2)]);
                                    clear tempBounds;
                                    
                                    %grey the current confidence text
                                    %set the confidence text
                                    for i=1:6
                                        Screen('drawtext',mainWindow,num2str(i),round(center_conf_text-(imageSize-conf_scale_bound)/2+((i-1)*((imageSize-conf_scale_bound)/5))-5),round(bottomY+imageSize/2+textImage_dist+textSpacing),textColorFade);
                                    end
                                    %
                                    tempBounds = Screen('TextBounds',mainWindow,confInst{4}); %followed?
                                    Screen('drawtext',mainWindow,confInst{4},round(center_conf_text-tempBounds(3)/2),round(bottomY+imageSize/2+textImage_dist-tempBounds(4)/5),textColorFade);
                                    clear tempBounds;
                                    
                                    tempBounds = Screen('TextBounds',mainWindow,confInst{1});
                                    Screen('drawtext',mainWindow,confInst{1},round(center_conf_text-(imageSize-conf_scale_bound)/2-tempBounds(3)/2),round(bottomY+imageSize/2+textImage_dist-tempBounds(4)/5+(2*textSpacing)),textColorFade);
                                    clear tempBounds;
                                    
                                    tempBounds = Screen('TextBounds',mainWindow,confInst{2});
                                    Screen('drawtext',mainWindow,confInst{2},round(center_conf_text-(imageSize-conf_scale_bound)/2-tempBounds(3)/2),round(bottomY+imageSize/2+textImage_dist-tempBounds(4)/5+(3*textSpacing-3)),textColorFade);
                                    clear tempBounds;
                                    
                                    tempBounds = Screen('TextBounds',mainWindow,confInst{1});
                                    Screen('drawtext',mainWindow,confInst{1},round(center_conf_text+(imageSize-conf_scale_bound)/2-tempBounds(3)/2),round(bottomY+imageSize/2+textImage_dist-tempBounds(4)/5+(2*textSpacing)),textColorFade);
                                    clear tempBounds;
                                    
                                    tempBounds = Screen('TextBounds',mainWindow,confInst{3});
                                    Screen('drawtext',mainWindow,confInst{3},round(center_conf_text+(imageSize-conf_scale_bound)/2-tempBounds(3)/2),round(bottomY+imageSize/2+textImage_dist-tempBounds(4)/5+(3*textSpacing-3)),textColorFade);
                                    clear tempBounds;
                                    
                                    %turn it al grey and update response time
                                    resp_start=Screen('Flip', mainWindow,0,1);
                                    
                                    %wait before continue - stuipd hack so that ptb
                                    %won't count the response of choosing the stimuli
                                    %as the response to the confidence
                                    while GetSecs<resp_start+0.2; end
                                  
                                    %get out of the confidence loop
                                    gave_conf_resp=1;
                                    
                                end
                                
                            end %ends if a conf response was given
                            
                        end %ends the while loop for duration of the trial
                    end %ends the conditional of wether the subject gave a response for the distractor
                    
                end
            end %ends if a distractor was chosen
            
            if (curr_trial_responses~=0) %subject gave 3 responses - move to the next trial
                removed=1;
                %clear the screen
                Screen(mainWindow,'FillOval',fixColor,fixDotRect);
                Screen('Flip',mainWindow);
            end
        end %ends the while loop for duration of the trial
        
        % fixation until next trial
        fix_start=GetSecs;
        while(GetSecs< fix_start + iti)
            WaitSecs(.0005);
            Screen(mainWindow,'FillOval',fixColor,fixDotRect);
            Screen('Flip',mainWindow);
        end
        
        % print output to screen and text file
        if practice==0
        fprintf('%d\t%s\t%d\t%d\t%.3f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.3f\t%.3f\t%.3f\t%d\t%d\t%d\t%.3f\t%.3f\t%.3f\n', ...
            subj_num, subj_id, trial, condition(trial), stim_onset(trial), ExpTest_all_trials_InLoc(trial,4), curr_trial_dist(trial,:),curr_trial_loc(trial,:),...
            conf_bt(trial,:),conf_rt(trial,:),choice_bt(trial,:), choice_rt(trial,:));
        
        fprintf(fid,'%d\t%s\t%d\t%d\t%.3f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.3f\t%.3f\t%.3f\t%d\t%d\t%d\t%.3f\t%.3f\t%.3f\n', ...
            subj_num, subj_id, trial, condition(trial), stim_onset(trial), ExpTest_all_trials_InLoc(trial,4), curr_trial_dist(trial,:),curr_trial_loc(trial,:),...
            conf_bt(trial,:),conf_rt(trial,:),choice_bt(trial,:), choice_rt(trial,:));
        end
        % plot performance and rt ongoing
%         if fMRI==1 || debug==1
%             rd.log_trial(response_times(run,trial), correct_response(run,trial), response(run,trial));
%         end
    end
    
    % wait for scanner to finish
    %while GetSecs < run_length; end
    
    % close this run's text file
    if practice==0
        fclose(fid);
    end
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

end % run loop

end

%% clean up and go home
if practice==0
    mat_file=fullfile(data_dir,'explicitPart2.mat');
    if exist(mat_file,'file')
       mat_file=[mat_file(1:end-4) '1.mat'];
    end
    save(mat_file,'condition','stim_onset','ExpTest_all_trials_InLoc','curr_trial_dist','curr_trial_loc','conf_bt','conf_rt','choice_bt','choice_rt','description')
end

if ~debug
    ListenChar(0);
    ShowCursor;
end
sca