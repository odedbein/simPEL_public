function [data_learning, simPEL9_All_subs_learning_byItemAllReps,simPEL9_All_subs_learning, simPEL9_All_subs_learning_averages,simPEL9_All_subs_learning_byItemAllRepsCollapsVnv]=simPEL9_analyse_AllLearningPhases_AB_mem_Oded()

project_dir='/Volumes/data/Bein/simPEL/simPEL9_onlyConsSameSimilarHalf';


analysis_dir=fullfile(project_dir,'analysis');
data_dir=fullfile(project_dir,'data');
addpath(genpath(data_dir));

%all subjs:
subjects_names={'CB','YP','SP','CB','AG','AS','MO','SS','KG','RK', 'PT', 'JT',...
                'YL', 'BE', 'OM', 'RJ', 'AW', 'TW', 'KB','GG','CL','TL','AB',...
                'CW','JM','NT','CJ','TW','PR','MG','JS','AA','EF','VL','NN','GH',...
                'JC','AK','IC'};
subjects_numbers=[2 4 5:8 12:18 20:30 34 36:39 41 43:48 80:82];
reps_init=9;


for subj=1:numel(subjects_numbers)
    subj_num=subjects_numbers(subj);
    subj_name=['S' num2str(subj_num) subjects_names{subj}];
    fprintf('currently analyzing subject %s \n',subj_name);
    data_learning.(subj_name)=analyze_encoding_single_sub(subj_num,subjects_names{subj},project_dir);
   
end

%produce tables with all subjects in each condition averages
part_name={'day1_init','rem','er_block'};
ABmem_name={'AllItems','Brem','Bforg'};
prev_chunk_size=2;
simPEL9_All_subs_learning={};
simPEL9_All_subs_learning_averages={};
simPEL9_All_subs_learning_byItemAllRepsCollapsVnv={};
simPEL9_All_subs_learning_byItemAllReps={};
for bb=1:numel(ABmem_name)
    for part=1:numel(part_name)

        if part==1
            rep_num=reps_init;
        else
            rep_num=1;
        end

        for rep=1:rep_num
            rep_name=['rep' num2str(rep)];
            for subj=1:length(subjects_numbers)
                subj_num=subjects_numbers(subj);
                subj_name=['S' num2str(subj_num) subjects_names{subj}];
                simPEL9_All_subs_learning.accuracy_num.(ABmem_name{bb}).(part_name{part}).(rep_name)(subj,:)=data_learning.(subj_name).(ABmem_name{bb}).(part_name{part}).accuracy_num(rep,:);
                simPEL9_All_subs_learning.accuracy_rates.(ABmem_name{bb}).(part_name{part}).(rep_name)(subj,:)=data_learning.(subj_name).(ABmem_name{bb}).(part_name{part}).accuracy_rates(rep,:);
                simPEL9_All_subs_learning.RT.(ABmem_name{bb}).(part_name{part}).(rep_name)(subj,:)=data_learning.(subj_name).(ABmem_name{bb}).(part_name{part}).RT(rep,:);
                simPEL9_All_subs_learning.logRT.(ABmem_name{bb}).(part_name{part}).(rep_name)(subj,:)=data_learning.(subj_name).(ABmem_name{bb}).(part_name{part}).logRT(rep,:);
            end
            simPEL9_All_subs_learning_averages.accuracy_num.(ABmem_name{bb}).(part_name{part})(rep,:)=mean(simPEL9_All_subs_learning.accuracy_num.(ABmem_name{bb}).(part_name{part}).(rep_name));
            simPEL9_All_subs_learning_averages.accuracy_rates.(ABmem_name{bb}).(part_name{part})(rep,:)=mean(simPEL9_All_subs_learning.accuracy_rates.(ABmem_name{bb}).(part_name{part}).(rep_name));
            simPEL9_All_subs_learning_averages.RT.(ABmem_name{bb}).(part_name{part})(rep,:)=mean(simPEL9_All_subs_learning.RT.(ABmem_name{bb}).(part_name{part}).(rep_name));
            simPEL9_All_subs_learning_averages.logRT.(ABmem_name{bb}).(part_name{part})(rep,:)=mean(simPEL9_All_subs_learning.logRT.(ABmem_name{bb}).(part_name{part}).(rep_name));
        end

    end
end

for subj=1:length(subjects_numbers)
    subj_num=subjects_numbers(subj);
    subj_name=['S' num2str(subj_num) subjects_names{subj}];
    
    %collapse on v and nv during initial learning and reminder
    for bb=1:numel(ABmem_name)
        for part=1:2

            if part==1
                rep_num=reps_init;
            else
                rep_num=1;
            end

            for loc=1:prev_chunk_size
                curr_trials=sum(data_learning.(subj_name).(ABmem_name{bb}).(part_name{part}).accuracy_num(:,[loc,loc+prev_chunk_size]),2);
                simPEL9_All_subs_learning_byItemAllRepsCollapsVnv.accuracy_num.(ABmem_name{bb}).(part_name{part})(subj,((loc-1)*rep_num)+(1:rep_num))=curr_trials';
                curr_trials=mean(data_learning.(subj_name).(ABmem_name{bb}).(part_name{part}).accuracy_rates(:,[loc,loc+prev_chunk_size]),2);
                simPEL9_All_subs_learning_byItemAllRepsCollapsVnv.accuracy_rates.(ABmem_name{bb}).(part_name{part})(subj,((loc-1)*rep_num)+(1:rep_num))=curr_trials';
                curr_trials=mean(data_learning.(subj_name).(ABmem_name{bb}).(part_name{part}).RT(:,[loc,loc+prev_chunk_size]),2);
                simPEL9_All_subs_learning_byItemAllRepsCollapsVnv.RT.(ABmem_name{bb}).(part_name{part})(subj,((loc-1)*rep_num)+(1:rep_num))=curr_trials';
                curr_trials=mean(data_learning.(subj_name).(ABmem_name{bb}).(part_name{part}).logRT(:,[loc,loc+prev_chunk_size]),2);
                simPEL9_All_subs_learning_byItemAllRepsCollapsVnv.logRT.(ABmem_name{bb}).(part_name{part})(subj,((loc-1)*rep_num)+(1:rep_num))=curr_trials';
            end

        end

        %this is the non-collapsed version
        for part=1:3

            curr_size=size(data_learning.(subj_name).(ABmem_name{bb}).(part_name{part}).accuracy_num);
            for item=1:curr_size(2) %for each item
                num_rep=curr_size(1);
                simPEL9_All_subs_learning_byItemAllReps.accuracy_num.(ABmem_name{bb}).(part_name{part})(subj,(item-1)*num_rep+(1:num_rep))=data_learning.(subj_name).(ABmem_name{bb}).(part_name{part}).accuracy_num(:,item)';
                simPEL9_All_subs_learning_byItemAllReps.accuracy_rate.(ABmem_name{bb}).(part_name{part})(subj,(item-1)*num_rep+(1:num_rep))=data_learning.(subj_name).(ABmem_name{bb}).(part_name{part}).accuracy_rates(:,item)';
                simPEL9_All_subs_learning_byItemAllReps.RT.(ABmem_name{bb}).(part_name{part})(subj,(item-1)*num_rep+(1:num_rep))=data_learning.(subj_name).(ABmem_name{bb}).(part_name{part}).RT(:,item)';
                simPEL9_All_subs_learning_byItemAllReps.logRT.(ABmem_name{bb}).(part_name{part})(subj,(item-1)*num_rep+(1:num_rep))=data_learning.(subj_name).(ABmem_name{bb}).(part_name{part}).logRT(:,item)';
            end

        end
    end
    
end

cd(fullfile(analysis_dir,'scripts'));
end

function curr_sub_data=analyze_encoding_single_sub(subj_num,subj_id,project_dir)
%IMPORTANT: this version had pairs instead of triads - I kept all the variable
%names including "triads" just because I was lazy to change them, but these
%refer to pairs.
prev_chunk_size=2; %in this version - pairs, in previous triads - would be 3.
data_dir=fullfile(project_dir,'data',sprintf('%d%s',subj_num,subj_id));
cd(data_dir);
part_name={'day1_init','rem','er_block'};
ABmem_name={'AllItems','Brem','Bforg'};
load([data_dir '/items_randomization.mat']);
load([data_dir '/trial_sequences.mat']);
%conditions:
%1-v 2-nv
%3-v violation 4-nv no violation (only in the error block)
%5-violation items
%6-no violation items

triads_per_cond=36; %number of pairs in each condition
loc_all=[1 2];
loc_er=1;

Exclude_SD=3; %exclusion criteria for RT outliers
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
    too_quick=(response_times<250);%find too quick responses (less then 250ms) - these are mistakes, remove them
    response_times(too_quick)=0;
    response_accuracy(too_quick)=0; %too quick responses shouldn't be counted as accurate.
    response_logRT=log(response_times);
    response_logRT(response_times==0)=0;
    
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
        curr_RTs(loc_in_pair(rep,:)==1)=0; %null A items - don't want to consider them because they are at the begining of a pair and the big/small responses are not controlled for.
        curr_RTs(response_accuracy(rep,:)==0)=0; %inaccurate responses - don't consider them as well, we don't know what the subjects did in them;
        curr_std=std(curr_RTs(curr_RTs~=0));
        curr_av=mean(curr_RTs(curr_RTs~=0));
        
        slow_outliers=find(curr_RTs>(curr_av+(Exclude_SD*curr_std)));
        fast_outliers=find(curr_RTs<(curr_av-(Exclude_SD*curr_std)));
        fast_outliers=fast_outliers(curr_RTs(fast_outliers)~=0);
        outliers=sort([slow_outliers fast_outliers]);
        if ~isempty(outliers)
            fprintf('outliers %s, rep %d: \n',part_name{part},rep);
            outliers
            sprintf('\n');
        end
        %mark these outliers:
        All_outliers.(part_name{part})(rep,outliers)=1;
    end
    
    
    %analyze based on ABmemory
    for bb = 1:numel(ABmem_name)
        curr_part=fullfile(data_dir,'explicitPart2.mat');
        load(curr_part,'ExpTest_all_trials_InLoc','choice_bt');
        curr_part=fullfile(data_dir,'items_randomization.mat');
        item_rand=load(curr_part);
        if bb == 1
           Bmem_in_enc=ones(size(response')); %take all trials - need to transpose so that it'll run on the columns
        else
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
           
        end
        %check that Bmem worked well:
        if any(isnan(Bmem_in_enc))
            fprintf('Bmem_in_enc has nans!!\n');
        end
        %transpose Bmem_in_enc to be the same dim as the response:
        Bmem_in_enc=Bmem_in_enc';
        %%%%% I STOPPED HERE
        if part<3
            accuracy_rates=zeros(size(response,1),num_triads_conds*prev_chunk_size);  %each condition * A B
        else
            accuracy_rates=zeros(size(response,1),num_triads_conds*prev_chunk_size+2); %each condition * A B + each condition single
        end
        accuracy_num=zeros(size(accuracy_rates));
        RT=zeros(size(accuracy_rates));
        logRT=zeros(size(accuracy_rates));
        
        %analyse v/nv: need to take into account the location in triad
        for cond=1:num_triads_conds
            c=parts_cond(cond);
            if (part==3 && c==3) %this is a violation condition in the vnv block - no 3 locs
                curr_loc=loc_er;
            else
                curr_loc=loc_all;
            end
            for loc=curr_loc %for all locs
                
                %select curr_trials based on loc, condition, and AB memory
                if bb==3 %forgoten AB
                   curr_trials=((loc_in_pair==loc)+(condition==c)+(Bmem_in_enc==0)==3);
                else
                   curr_trials=((loc_in_pair==loc)+(condition==c)+(Bmem_in_enc==1)==3); %remembered, or, if take all, all will be marked as "1"
                end
                temp=zeros(size(response));
                temp_all_resp=zeros(size(response));
                temp(curr_trials)=response_accuracy(curr_trials);
                accuracy_rates(:,(cond-1)*prev_chunk_size+loc)=sum(temp,2)./sum(curr_trials,2);
                accuracy_num(:,(cond-1)*prev_chunk_size+loc)=sum(temp,2);
                temp=zeros(size(response));
                temp(curr_trials)=response_times(curr_trials);
                temp(All_outliers.(part_name{part})==1)=0; %null if outlier - if not the current condition, anyway it's zero, no problem.
                temp_all_resp(curr_trials)=response(curr_trials)~=0;
                %to know how many trials were there eventually, after outliers
                %removal and removal of inaccurate responses:
                temp_num_no_outliers=temp;
                temp_num_no_outliers(temp>0)=1;
                temp_num_no_outliers(response_accuracy==0)=0;
                temp_num_no_outliers=sum(temp_num_no_outliers,2);
                
                %if subj didn't reply he'll have 0 in his RT mat, so no need to
                %remove it, but if they replied but inaccurate - I'll have it.
                %so do remove - but not for the first location
                if loc == 1 %This is an A item, don't look at accuracy, it's meaningless, only remove responses if subject did not give a response at all
                    num_responses=sum(temp_all_resp,2); %I don't remove outliers from the A items, so can take this rather than temp_num_no_outliers
                    RT(:,(cond-1)*prev_chunk_size+loc)=sum(temp,2)./num_responses;
                else
                    temp(response_accuracy==0)=0; %only for B items,consider accuracy A item RT is calculated regardless of accuracy, because it is meaningless
                    RT(:,(cond-1)*prev_chunk_size+loc)=sum(temp,2)./temp_num_no_outliers;
                end
                
                temp=zeros(size(response));
                temp(curr_trials)=response_logRT(curr_trials);
                temp(All_outliers.(part_name{part})==1)=0; %null if outlier - if not the current condition, anyway it's zero, no problem
                % NOTE: outliers are computed based on the data in ms - i think
                % it makes sense, because this is the measure in which the
                % participant response. log RT is just for stats etc.
                
                if loc == 1 %This is an A item, don't look at accuracy, it's meaningless, only remove responses if subject did not give a response at all
                    num_responses=sum(temp_all_resp,2);
                    logRT(:,(cond-1)*prev_chunk_size+loc)=sum(temp,2)./num_responses;
                else
                    temp(response_accuracy==0)=0; %only for B and C items,consider accuracy A item RT is calculated regardless of accuracy, because it is meaningless
                    logRT(:,(cond-1)*prev_chunk_size+loc)=sum(temp,2)./temp_num_no_outliers;
                end
            end
        end
        
        %analyse single items, only vnv block, and no location in triad
        if part == 3
            for cond=num_triads_conds+1:length(parts_cond)
                c=parts_cond(cond);
                if bb==3 %forgoten AB
                   curr_trials=((condition==c)+(Bmem_in_enc==0)==2);
                else
                   curr_trials=((condition==c)+(Bmem_in_enc==1)==2); %remembered, or, if take all, all will be marked as "1"
                end
                temp=zeros(size(response));
                temp(curr_trials)=response_accuracy(curr_trials);
                accuracy_rates(:,num_triads_conds*prev_chunk_size+c-num_triads_conds)=sum(temp,2)./sum(curr_trials,2);
                accuracy_num(:,num_triads_conds*prev_chunk_size+c-num_triads_conds)=sum(temp,2);
                temp=zeros(size(response));
                temp(curr_trials)=response_times(curr_trials);
                temp(All_outliers.(part_name{part})==1)=0; %null if outlier - if not the current condition, anyway it's zero, no problem
                %to know how many trials were there eventually, after outliers
                %removal:
                temp_num_no_outliers=temp;
                temp_num_no_outliers(temp>0)=1;
                temp_num_no_outliers(response_accuracy==0)=0;
                temp_num_no_outliers=sum(temp_num_no_outliers,2);
                
                %if subj didn't reply I'll have 0 in his RT mat, so no need to
                %remove it, but if they replied but inaccurate - I'll have it.
                %so do remove
                temp(response_accuracy==0)=0;
                RT(:,num_triads_conds*prev_chunk_size+c-num_triads_conds)=sum(temp,2)./temp_num_no_outliers;
                temp=zeros(size(response));
                temp(curr_trials)=response_logRT(curr_trials);
                temp(response_accuracy==0)=0;
                temp(All_outliers.(part_name{part})==1)=0; %null if outlier - if not the current condition, anyway it's zero, no problem
                logRT(:,num_triads_conds*prev_chunk_size+c-num_triads_conds)=sum(temp,2)./temp_num_no_outliers;
            end
        end
        
        %save it all:
        curr_sub_data.(ABmem_name{bb}).(part_name{part}).accuracy_num=accuracy_num;
        curr_sub_data.(ABmem_name{bb}).(part_name{part}).accuracy_rates=accuracy_rates;
        curr_sub_data.(ABmem_name{bb}).(part_name{part}).RT=RT;
        curr_sub_data.(ABmem_name{bb}).(part_name{part}).logRT=logRT;
        
    end %ends the ABmem
    
end

end

