function data_learning=simPEL9_compile_encodingRT(engram) 
if engram
    mydir='/data/Bein';
else
    mydir='/Volumes/data/Bein';
end

project_dir=fullfile(mydir,'simPEL/simPEL9_onlyConsSameSimilarHalf');
analysis_dir=fullfile(project_dir,'analysis');
results_dir=fullfile(analysis_dir,'files_forR');
if ~exist(results_dir); mkdir(results_dir); end
data_dir=fullfile(project_dir,'data');
addpath(genpath(data_dir));

%all subjs:
subjects_names={'CB','YP','SP','CB','AG','AS','MO','SS','KG','RK', 'PT', 'JT',...
                'YL', 'BE', 'OM', 'RJ', 'AW', 'TW', 'KB','GG','CL','TL','AB',...
                'CW','JM','NT','CJ','TW','PR','MG','JS','AA','EF','VL','NN','GH',...
                'JC','AK','IC'};
subjects_numbers=[2 4 5:8 12:18 20:30 34 36:39 41 43:48 80:82];

good_subj=[2 3 6 8:17 19:20 22 24:26 31:39];
%to make things easier, let's exclude them from now:
subjects_names=subjects_names(good_subj);
subjects_numbers=subjects_numbers(good_subj);


for subj=1:numel(subjects_numbers)
    subj_num=subjects_numbers(subj);
    subj_name=['S' num2str(subj_num) subjects_names{subj}];
    fprintf('currently analyzing subject %s \n',subj_name);
    data_learning.(subj_name)=analyze_encoding_single_sub(subj_num,subjects_names{subj},project_dir);
end

%produce tables with all subjects in each condition averages
part_name={'day1_init','rem','er_block'};
header={'subject' 'repetition' 'condition' 'loc_in_pair' 'resp_big_small' 'resp_acc' 'RTs' 'Bmemory'};
for part=1:numel(part_name)
    all_sub_mat=[];
    for subj=1:length(subjects_numbers)
        subj_num=subjects_numbers(subj);
        subj_name=['S' num2str(subj_num) subjects_names{subj}];
        all_sub_mat=[all_sub_mat;data_learning.(subj_name).(part_name{part})];
    end
    
    % set the table:
    T=array2table(all_sub_mat);
    T.Properties.VariableNames=header;
    %write it up:
    results_fname=['simPEL9_associative_encoding_' part_name{part} '.txt'];
    filename=fullfile(results_dir,results_fname);
    writetable(T,filename,'Delimiter','\t')

end
cd(fullfile(analysis_dir,'scripts'));
end

function curr_sub_data=analyze_encoding_single_sub(subj_num,subj_id,project_dir)

data_dir=fullfile(project_dir,'data',sprintf('%d%s',subj_num,subj_id));
cd(data_dir);
part_name={'day1_init','rem','er_block'};
load([data_dir '/items_randomization.mat']);
load([data_dir '/trial_sequences.mat']);
%conditions:
%1-v 2-nv
%3-v violation 4-nv no violation (only in the error block)
%5-violation items
%6-no violation items

%for the RT - criteria for exclusion:
Exclude_SD=3;
too_quick_thresh=250;

All_outliers={}; %to save RT outliers

for part=1:numel(part_name)
    curr_part=[part_name{part} '_' num2str(subj_num) subj_id '.mat'];
    load(curr_part,'response','response_times','response_accuracy','condition','stim_in_pair');
    
    if part<3 %meaning - not the vnv block
        parts_cond=[1 2]; %which conditions are in this part - see condition list above
        num_triads_conds=2; %how many conditions have pairs - important for later.
    else
        parts_cond=1:6; %which conditions are in this part - see condition list above
        num_triads_conds=4;%how many conditions have pairs. conditions 5,6 are the v/nv items, so they are not pairs
    end
    
    switch part
        case 1
            curr_seq=stim_seq_day1;
            loc_in_pair=stim_loc_seq_day1;
            stim_in_pair=stim_in_pair_seq_day1';
        case 2
            curr_seq=stim_seq_reminder;
            loc_in_pair=stim_loc_seq_reminder;
            stim_in_pair=stim_in_pair_seq_reminder';
        case 3
            curr_seq=stim_seq_vnv;
            loc_in_pair=stim_loc_seq_vnv;
            stim_in_pair=stim_in_pair_seq_vnv';
    end
    
    %prep matrices
    
    response_times=response_times*1000;
    too_quick=(response_times<too_quick_thresh);%find too quick responses (less then 250ms) - these are mistakes, remove them
    response_times(too_quick)=nan;
    response_accuracy(isnan(too_quick))=0; %too quick responses shouldn't be counted as accurate.
    response_logRT=log(response_times);
    response_logRT(isnan(response_times))=0;
    
    %in the initial learning phase, each run is a repetition, and outliers are computed per run (later). In the reminder and the violation block, we just gave ss breaks, runs are not repetitions, so
    %treat them all as one repetition.
    if part > 1
        response=reshape(response',1,numel(response));
        response_accuracy=reshape(response_accuracy',1,numel(response_accuracy));
        response_times=reshape(response_times',1,numel(response_times));
        response_logRT=reshape(response_logRT',1,numel(response_logRT));
        loc_in_pair=reshape(loc_in_pair',1,numel(loc_in_pair));
        curr_seq=reshape(curr_seq',1,numel(curr_seq'));%I later on use it based on loc_in_pair, so if I rearrange, I need to rearrange this as well
        condition=reshape(condition',1,numel(condition));
        
    end
    
    %exclude RT outliers - each repetition has enough items, so do for each rep (previous versions I did on each pair of repetitions):
    All_outliers.(part_name{part})=zeros(size(response));
    
    for rep=1:size(response,1)
        curr_RTs=response_times(rep,:);
        curr_RTs(loc_in_pair(rep,:)==1)=nan; %null A items - don't want to consider them because they are at the begining of a pair and the big/small responses are not controlled for.
        curr_RTs(response_accuracy(rep,:)==0)=nan; %inaccurate responses - don't consider them as well, we don't know what the subjects did in them;
        curr_std=nanstd(curr_RTs);
        curr_av=nanmean(curr_RTs);
        
        slow_outliers=find(curr_RTs>(curr_av+(Exclude_SD*curr_std)));
        fast_outliers=find(curr_RTs<(curr_av-(Exclude_SD*curr_std)));
        fast_outliers=fast_outliers(curr_RTs(fast_outliers)~=0);
        outliers=sort([slow_outliers fast_outliers]);
        if ~isempty(outliers)
            fprintf('outliers %s, rep %d: \n',part_name{part},rep);
            outliers
            fprintf('\n')
        end
        %mark these outliers:
        All_outliers.(part_name{part})(rep,outliers)=1;
    end
    
    
    %analyze ABmemory
    
    curr_part=fullfile(data_dir,'explicitPart2.mat');
    load(curr_part,'ExpTest_all_trials_InLoc','choice_bt');
    curr_part=fullfile(data_dir,'items_randomization.mat');
    item_rand=load(curr_part);
    
    Bmem_in_enc=nan(size(response')); %initiate with zeros, to be filled where AB is rememebered.
    %since I sometime reshape loc
    Aitem_enc=unique(curr_seq(loc_in_pair==1));
    
    %Aitem_enc=[Aitem_enc nan(size(Aitem_enc))];
    for ai=1:numel(Aitem_enc)
        curr_cue=Aitem_enc(ai);
        all_Aloc=find(curr_seq'==curr_cue);
        %in the violation part, Aitem_enc was sometimes starting with
        %14/13, to denote the condition.
        %need to modify that to find in the memory test.
        if part==3
            if ismember(floor(abs(curr_cue/1000)),[13,14,23,24])
                if curr_cue < 0
                    curr_cue=curr_cue+2000;
                else
                    curr_cue=curr_cue-2000;
                end
            end
        end
        Bresp=choice_bt((ExpTest_all_trials_InLoc(:,4)==curr_cue),1)>0;
        if isempty(Bresp)
            fprintf('didn''t find Bresp!!\n');
        else
            Bmem_in_enc(all_Aloc)=Bresp;%put the response in the A loc
            Bmem_in_enc(all_Aloc+1)=Bresp;%put the response in the B loc
            if part == 3  %violation pairs - need to put that for the no-violation
                more_loc=all_Aloc+2;
                if more_loc <= length(Bmem_in_enc) %it's only in part 3, Bmem is a vector
                    add_loc=find(loc_in_pair(more_loc)==4);
                    Bmem_in_enc(more_loc(add_loc))=Bresp;
                end
            end
        end
    end
    
    %check that Bmem worked well:
    if any(isnan(Bmem_in_enc))
        fprintf('Bmem_in_enc has nans!!\n');
    end
    %transpose Bmem_in_enc to be the same dim as the response:
    Bmem_in_enc=Bmem_in_enc';
    
    %% build the matrix:
    %in some parts of the task the reshape/transpose might be redundant,
    %but it doesn't hurt, so for ease, I left it there.
    %condition (v/nv):
    condition_col=condition';
    condition_col=reshape(condition_col,numel(condition_col),1);
    %loc (first/second item in pair):
    loc_in_pair_col=loc_in_pair';
    loc_in_pair_col=reshape(loc_in_pair_col,numel(loc_in_pair_col),1);
    %repetition:
    repetition=repmat(1:size(condition,1),size(condition,2),1);
    repetition_col=reshape(repetition,numel(repetition),1);
    
    %response (bigger/smaller):
    resp_bg=nan(size(response));
    if (mod(subj_num,2) == 0) %even number - J is bigger, K is smaller. The numbers might differ based
        %on the computer that the participants ran, but it's J is always a
        %smaller number than K. I put here the exact numbers because if a
        %subject pressed h/l - I don't want to include it.
        if any(any(response==74))
            bigger_key=74;
            smaller_key=75;
        else
           bigger_key=13;
           smaller_key=14;
        end
    else
        if any(any(response==74))
            bigger_key=75;
            smaller_key=74;
        else
           bigger_key=14;
           smaller_key=13;
        end
    end
    
    resp_bg(response==bigger_key)=1;
    resp_bg(response==smaller_key)=0;
    resp_bg_col=reshape(resp_bg',numel(resp_bg),1);
    
    if all(isnan(resp_bg_col)==1)
        fprintf('subj %d%s part %d: smt is wrong, nan in bg responses\n',subj_num, subj_id,part)
    end
    
    %response accuracy (irrelevant for first items):
    resp_acc_col=reshape(response_accuracy',numel(response_accuracy),1);
    %RTs (in ms, nan outliers)
    RTs=response_times;
    RTs(response==0)=nan; %nan no responses
    RTs(All_outliers.(part_name{part})==1)=nan;
    RTs_col=reshape(RTs',numel(RTs),1);
    %Bmemory:
    Bmem_col=reshape(Bmem_in_enc',numel(Bmem_in_enc),1);
    subj_col=repmat(subj_num,numel(Bmem_in_enc),1);
    %save it all in one matrix:
    subj_mat=[subj_col repetition_col condition_col loc_in_pair_col resp_bg_col resp_acc_col RTs_col Bmem_col];
    curr_sub_data.(part_name{part})=subj_mat;
    
end



end

