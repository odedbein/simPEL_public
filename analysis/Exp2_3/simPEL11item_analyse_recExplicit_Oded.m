function [simPEL11_All_subs_recognition, simPEL11_All_subs_explicit]=simPEL11item_analyse_recExplicit_Oded()

%%
project_dir='/Volumes/data/Bein/simPEL/simPEL11_onlyConsSameSimilarItemTask';

exclude_outliers=0;%THERE's A BUG: this because of the loop of the runs, it also influences accuracy data = 
%check if you want to remove outliers. I didn't in all analyses reported in the paper, so doesn't matter. only from RT data


%all subjs:
subjects_names={'JS','JP','KE','PK','EF','CM','MH','OH','JD','SM','AM','CM','SG','JP','EF','SE','SR','JB','NR','KZ','SC','RK','RY','NL','EK','BD','LT','MI'};
subjects_numbers=[1:15 18:30];

total_subj_num=numel(subjects_numbers);

%set up a figure for marking outliers:
f = figure;
set(f,'name','simPEL11','numbertitle','off');

%set up data structures:
simPEL11_All_subs_recognition={};
simPEL11_All_subs_explicit={};


for subj=1:numel(subjects_numbers)
    curr_sub=subjects_numbers(subj);
    subj_name=['S' num2str(curr_sub) subjects_names{subj}];
    fprintf('currently analyzing subject %s \n',subj_name);
    
    %analyze recognition
    curr_sub_data=analyze_recognition_single_sub(curr_sub,subjects_names{subj},project_dir,subj,total_subj_num,exclude_outliers);

    simPEL11_All_subs_recognition.All.num_items(subj,:)=curr_sub_data.All.num_items;
    simPEL11_All_subs_recognition.All.num_itemsRT(subj,:)=curr_sub_data.All.num_itemsRT;
    simPEL11_All_subs_recognition.All.total_responses(subj,:)=curr_sub_data.All.total_responses;
    simPEL11_All_subs_recognition.All.accuracy_num(subj,:)=curr_sub_data.All.accuracy_num;
    simPEL11_All_subs_recognition.All.accuracy_rates(subj,:)=curr_sub_data.All.accuracy_rates;
    simPEL11_All_subs_recognition.All.response_typePerCond(subj,:)=curr_sub_data.All.response_typePerCond;
    simPEL11_All_subs_recognition.All.response_OldOrSim(subj,:)=curr_sub_data.All.response_OldOrSim;
    simPEL11_All_subs_recognition.All.RT(subj,:)=curr_sub_data.All.RT;
    simPEL11_All_subs_recognition.All.logRT(subj,:)=curr_sub_data.All.logRT;

    simPEL11_All_subs_recognition.All.Bmem.num_items(subj,:)=curr_sub_data.All.Bmem.num_items;
    simPEL11_All_subs_recognition.All.Bmem.response_typePerCondNum(subj,:)=curr_sub_data.All.Bmem.response_typePerCondNum;
    simPEL11_All_subs_recognition.All.Bmem.response_typePerCond(subj,:)=curr_sub_data.All.Bmem.response_typePerCond;
     
    %analyze explicit
    ExpPart=1;
    [curr_sub_conf, curr_sub_confSD, curr_sub_OldNew, curr_sub_CountEachRating, ...
     curr_sub_confAvDist, curr_sub_confSDAvDist, curr_sub_OldNewAvDist, curr_sub_CountEachRatingAvDist]=analyze_explicit_single_sub(curr_sub,subjects_names{subj},project_dir,ExpPart);
    simPEL11_All_subs_explicit.Part1.AvConfRates(subj,:)=curr_sub_conf;
    simPEL11_All_subs_explicit.Part1.AvConfRatesSD(subj,:)=curr_sub_confSD;
    simPEL11_All_subs_explicit.Part1.OldNew(subj,:)=curr_sub_OldNew;
    simPEL11_All_subs_explicit.Part1.CountEachRating(subj,:)=curr_sub_CountEachRating;
    
    simPEL11_All_subs_explicit.Part1.AvConfRatesAvDist(subj,:)=curr_sub_confAvDist;
    simPEL11_All_subs_explicit.Part1.AvConfRatesSDAvDist(subj,:)=curr_sub_confSDAvDist;
    simPEL11_All_subs_explicit.Part1.OldNewAvDist(subj,:)=curr_sub_OldNewAvDist;
    simPEL11_All_subs_explicit.Part1.CountEachRatingAvDist(subj,:)=curr_sub_CountEachRatingAvDist;
    ExpPart=2;
    [curr_sub_conf, curr_sub_confSD, curr_sub_OldNew, curr_sub_CountEachRating, ...
     curr_sub_confAvDist, curr_sub_confSDAvDist, curr_sub_OldNewAvDist, curr_sub_CountEachRatingAvDist]=analyze_explicit_single_sub(curr_sub,subjects_names{subj},project_dir,ExpPart);
    simPEL11_All_subs_explicit.Part2.AvConfRates(subj,:)=curr_sub_conf;
    simPEL11_All_subs_explicit.Part2.AvConfRatesSD(subj,:)=curr_sub_confSD;
    simPEL11_All_subs_explicit.Part2.OldNew(subj,:)=curr_sub_OldNew;
    simPEL11_All_subs_explicit.Part2.CountEachRating(subj,:)=curr_sub_CountEachRating;
    
    simPEL11_All_subs_explicit.Part2.AvConfRatesAvDist(subj,:)=curr_sub_confAvDist;
    simPEL11_All_subs_explicit.Part2.AvConfRatesSDAvDist(subj,:)=curr_sub_confSDAvDist;
    simPEL11_All_subs_explicit.Part2.OldNewAvDist(subj,:)=curr_sub_OldNewAvDist;
    simPEL11_All_subs_explicit.Part2.CountEachRatingAvDist(subj,:)=curr_sub_CountEachRatingAvDist;
end
 
end


function curr_sub_data=analyze_recognition_single_sub(subj_num,subj_id,project_dir,subj,total_subj_num,exclude_outliers)
data_dir=fullfile(project_dir,'data',sprintf('%d%s',subj_num,subj_id));
curr_part=fullfile(data_dir,'recognition.mat');
load(curr_part,'response','response_times','response_accuracy','condition','priming_sequence','response_type');

response_times=response_times*1000;
nRuns=size(response_times,1);

%conditions:
%5-violation-old
%6-no-violation-old
%7-violation-similar
%8-no-violation-similar
%9-new

%response type:
%1-new
%2-similar
%3-old
%% analyze REC 
%header.rec1All={'violation','no violation','foils'};
parts_cond=[5 6 7 8 9];
num_cond_all=length(parts_cond);

too_quick=(response_times<250);%find too quick responses - these are mistakes, remove them
response_times(too_quick)=0;
response_accuracy(too_quick)=0;%if by chance they were accurate - remove them from accuracy.
response_type(too_quick)=0;

response_logRT=log(response_times);
response_logRT(response_times==0)=0;

%preparing the data structures:
accuracy_rates=zeros(nRuns,num_cond_all);
accuracy_num=zeros(nRuns,num_cond_all);
total_responses=zeros(nRuns,num_cond_all);
response_typePerCond=zeros(nRuns,num_cond_all*3); %for all conditions (5 conditions), all response types (1-3)
response_OldOrSim=zeros(nRuns,num_cond_all);
RT=zeros(nRuns,num_cond_all*3);
num_items=zeros(nRuns,num_cond_all);
num_itemsRT=zeros(nRuns,num_cond_all);
logRT=zeros(nRuns,num_cond_all*3);

%% mark outliers:
Exclude_SD=3; %exclusion criteria
RT_std=std(response_times(response_accuracy~=0));
RT_av=mean(response_times(response_accuracy~=0));
slow_outliers=find(response_times>(RT_av+(Exclude_SD*RT_std)));
slow_outliers=slow_outliers(response_accuracy(slow_outliers)~=0);%remove inaccurate responses from being marked as outliers. This does not matter for the analysis in which they'll be removed anyway, but more for mee to know how many real outleirs were there.
fast_outliers=find(response_times<(RT_av-(Exclude_SD*RT_std)));
fast_outliers=fast_outliers(response_accuracy(fast_outliers)~=0); %remove inaccurate responses from being marked as outliers. This does not matter for the analysis in which they'll be removed anyway, but more for mee to know how many real outleirs were there.
outliers=sort([slow_outliers;fast_outliers]);

%this part does not influence the analysis, only for plotting, so I'd plot only accurate responses:
only_accRT=response_times(response_accuracy~=0);
only_accCondition=condition(response_accuracy~=0);
outliers_new_loc=[];
for out=1:length(outliers)
    outliers_new_loc=[outliers_new_loc find(only_accRT==response_times(outliers(out)))];
end

%scatter plot:
subj_name=[num2str(subj_num) subj_id];
c=zeros(length(only_accRT),3);
c(outliers_new_loc,2)=1;
for out=outliers_new_loc
    if ismember(only_accCondition(out),[5 6]) %outlier is one of the critical conditions - change to red
        c(out,2)=0;
        c(out,1)=1;
    end
end
subplot(ceil(sqrt(total_subj_num)),ceil(sqrt(total_subj_num)), subj);
scatter(1:length(only_accRT),only_accRT,[],c)
ylabel('RTs (s)');
title(subj_name);

if ~isempty(outliers)
    outliers
    fprintf('\n');
    fprintf('conditions: %d \n',condition(outliers));
end


%now analyse the data:
for cond=1:num_cond_all
    c=parts_cond(cond);
    curr_trials=(condition==c); %marks the trials in the current condition
    num_items(:,cond)=sum(curr_trials,2);
    
    for r=1:nRuns
        accuracy_rates(r,cond)=sum(response_accuracy(r,curr_trials(r,:)))/num_items(r,cond);
        accuracy_num(r,cond)=sum(response_accuracy(r,curr_trials(r,:)));
        total_responses(r,cond)=sum(response(r,curr_trials(r,:))~=0);
        for t=1:3 %three response types
            response_typePerCond(r,(cond-1)*3+t)=sum(response_type(r,curr_trials(r,:))==t);
            
            %now deal with RTs - if want to exclude RT outliers from the
            %accuracy analysis - correct this! there's a bug! it'll remove
            %outliers only after calculating accuracy for the first run, first trial type. 
            if exclude_outliers
                curr_trials(outliers)=0; %remove RT outliers from the analysis, will do it for all runs, it's fine
            end
            curr_trials_temp=(curr_trials(r,:)+(response_type(r,:)==t))==2; %take only curr condition and also accurate responses
            
            temp=response_times(r,curr_trials_temp);
            RT(r,(cond-1)*3+t)=nansum(temp);
            num_itemsRT(r,(cond-1)*3+t)=length(temp);
            temp=response_logRT(curr_trials_temp);
            logRT(r,(cond-1)*3+t)=nansum(temp);
        end
        %OldOrSim: summarize across both old and sim responses: for each
        %participants, how many times did she responded either sim or old.
        response_OldOrSim(r,cond)=nansum((response_type(r,curr_trials(r,:))==2)+(response_type(r,curr_trials(r,:))==3));
        
        
    end
end

%save it all:
curr_sub_data.All.num_items=sum(num_items);
curr_sub_data.All.num_itemsRT=sum(num_itemsRT);
curr_sub_data.All.total_responses=sum(total_responses);
curr_sub_data.All.accuracy_num=sum(accuracy_num);
curr_sub_data.All.accuracy_rates=mean(accuracy_rates);
curr_sub_data.All.response_typePerCond=sum(response_typePerCond)./(reshape(repmat(sum(num_items),3,1),1,numel(repmat(sum(num_items),1,3))));
curr_sub_data.All.response_OldOrSim=sum(response_OldOrSim)./sum(num_items);
curr_sub_data.All.RT=nansum(RT)./nansum(num_itemsRT);
curr_sub_data.All.logRT=nansum(logRT)./nansum(num_itemsRT);

%% separate by memory for the B item 
parts_cond=[5 6 7 8];
num_cond_all=length(parts_cond);
response_typePerCond=zeros(nRuns,num_cond_all*3*2); %like previously, but w/o new items and by rem/for B items
total_responses=zeros(1,num_cond_all*2);%v/nv by old/similar by rem/for B items


curr_part=fullfile(data_dir,'explicitPart2.mat');
load(curr_part,'ExpTest_all_trials_InLoc','choice_bt');
curr_part=fullfile(data_dir,'items_randomization.mat');
load(curr_part);
curr_part=fullfile(data_dir,'trial_sequences.mat');
load(curr_part,'stim_seq_vnv');
viol_seq=stim_seq_vnv';
Bmem=zeros(size(priming_sequence));

%prep the Bmem matrix - says which items were later B remembered and which weren't 
for cond=1:2%length(parts_cond)
    %that's for the old items
    c=parts_cond(cond);
    curr_trials=(condition==c);
    curr_items=priming_sequence(curr_trials);
    for i=1:length(curr_items)
        item=curr_items(i);
        if cond < 2
            cue=viol_seq(find(viol_seq==item)-1);%a item was immediately before
        else
            cue=viol_seq(find(viol_seq==item)-2);%a item was 2 before before
        end
        if cue < 0 
            cue=cue+2000;
        else
            cue=cue-2000;
        end
        
        if isempty(find(ExpTest_all_trials_InLoc(:,4)==cue))
            display(sprintf('something is wrong, didn''t find the cue of item number %d, condition %d\n',i,c))
        else
            Bmem(priming_sequence==item)=(choice_bt((ExpTest_all_trials_InLoc(:,4)==cue),1)>0); %put the accuracy in cMem, subject chose the correct response
        end
    end
    
    %now for similar lures items
    c=parts_cond(cond+2);
    curr_trials=(condition==c);
    curr_items=priming_sequence(curr_trials); 
    for i=1:length(curr_items)
        %lures and old were always pairs of 2 stimuli - for each pair 1-2
        %either 1 was old and 2 was similar or the opposite. That means
        %that if it's an even number - the old one is one before, while if
        %it's an odd number - the old one is one after.
        if mod(curr_items(i),2)==0 %even
            item=curr_items(i)-2000-1;
        else %odd
            item=curr_items(i)-2000+1;
        end
        if cond < 2
            cue=viol_seq(find(abs(viol_seq)==item)-1);%a item was immediately before %need absolute values of the viol_seq because my item is positive
        else
            cue=viol_seq(find(abs(viol_seq)==item)-2);%a item was 2 before before
        end
        if cue < 0 
            cue=cue+2000;
        else
            cue=cue-2000;
        end
        
        if isempty(find(ExpTest_all_trials_InLoc(:,4)==cue))
            display(sprintf('something is wrong, didn''t find the cue of item number %d, condition %d\n',i,c))
        else
            Bmem(priming_sequence==curr_items(i))=(choice_bt((ExpTest_all_trials_InLoc(:,4)==cue),1)>0); %put the accuracy in cMem, subject chose the correct response
        end    
    end
end
    
for cond=1:length(parts_cond)    
    %remembered B items
    c=parts_cond(cond);
    curr_trials=(condition==c);
    curr_trials(Bmem==0)=0; %null all forgotten
    for r=1:nRuns
        total_responses(r,cond)=sum(response(r,curr_trials(r,:))~=0);
        for t=1:3 %three response types
            response_typePerCond(r,(cond-1)*3+t)=sum(response_type(r,curr_trials(r,:))==t);
        end
       
    end   
    
    %forgotten B items
    curr_trials=(condition==c);
    curr_trials(Bmem==1)=0; %null all remembered
    for r=1:nRuns
        total_responses(r,num_cond_all+cond)=sum(response(r,curr_trials(r,:))~=0);
        for t=1:3 %three response types
            response_typePerCond(r,num_cond_all*3+(cond-1)*3+t)=sum(response_type(r,curr_trials(r,:))==t);
        end
       
    end   
end

curr_sub_data.All.Bmem.num_items=sum(total_responses(1:nRuns,:));
curr_sub_data.All.Bmem.response_typePerCondNum=sum(response_typePerCond(1:nRuns,:));
curr_sub_data.All.Bmem.response_typePerCond=sum(response_typePerCond(1:nRuns,:))./...
                                               (reshape(repmat(curr_sub_data.All.Bmem.num_items,3,1),1,numel(repmat(curr_sub_data.All.Bmem.num_items,1,3))));
                                          
end

%% analyse explicit memory test
function [curr_sub_conf, curr_sub_confSD, curr_sub_Old, curr_sub_CountEachRating, curr_sub_confAvDist, curr_sub_confSDAvDist, curr_sub_OldAvDist, curr_sub_CountEachRatingAvDist] = analyze_explicit_single_sub(subj_num,subj_id,project_dir,ExpPart)
data_dir=fullfile(project_dir,'data',sprintf('%d%s',subj_num,subj_id));
curr_part=fullfile(data_dir,sprintf('explicitPart%d.mat',ExpPart));
load(curr_part);
num_ratings=6;%number of points on confidence scale
triads_per_cond=36;
%calculate average confidence rates of the subject: 
av_conf=reshape(conf_bt,[numel(conf_bt),1]);
av_conf=mean(av_conf(av_conf~=0)); %remove zeors - trials in which the subject did not respond
sd_conf=reshape(conf_bt,[numel(conf_bt),1]);
sd_conf=std(sd_conf(sd_conf~=0));%remove zeors - trials in which the subject did not respond
conf_btSTD=conf_bt-av_conf;
conf_btSTD=conf_btSTD/sd_conf;
conf_btSTD(conf_bt==0)=0;
%conditions:
%1-v 2-nv
condition=rem(ExpTest_all_trials_InLoc(:,4),10000);
condition=floor(abs(condition/1000));
num_cond=length(unique(condition));
curr_sub_conf=zeros(1,3*num_cond);
curr_sub_confSD=zeros(1,3*num_cond);
curr_sub_Old=zeros(1,3*num_cond);
curr_sub_CountEachRating=zeros(1,3*num_ratings*num_cond);

curr_sub_confAvDist=zeros(1,2*num_cond);
curr_sub_confSDAvDist=zeros(1,2*num_cond);
curr_sub_OldAvDist=zeros(1,2*num_cond);
curr_sub_CountEachRatingAvDist=zeros(1,2*num_ratings*num_cond);

for cond=1:num_cond
    curr_trials=(condition==cond);
    for loc=1:3
        curr_conf=conf_bt(curr_trials,loc);
        curr_conf=curr_conf(curr_conf~=0); %take only those in which the participant chose this answer (loc=1 is the target)
        curr_sub_conf((cond-1)*3+loc)=sum(curr_conf)/length(curr_conf); %average confidence rates`
        curr_sub_Old((cond-1)*3+loc)=length(curr_conf)/triads_per_cond; %accuracy rates
        for rate=1:num_ratings
            curr_sub_CountEachRating((cond-1)*3*num_ratings+(loc-1)*num_ratings+rate)=sum(curr_conf==rate);
        end
        curr_conf=conf_btSTD(curr_trials,loc);
        curr_conf=curr_conf(curr_conf~=0); %remove if subject did not respond
        curr_sub_confSD((cond-1)*3+loc)=sum(curr_conf)/length(curr_conf);
    end
end

%average all distractors
for cond=1:num_cond
    curr_trials=(condition==cond);
    loc=1;
        curr_conf=conf_bt(curr_trials,loc);
        curr_conf=curr_conf(curr_conf~=0); %take only those in which the participant chose this answer (loc=1 is the target)
        curr_sub_confAvDist((cond-1)*2+loc)=sum(curr_conf)/length(curr_conf); %average
        curr_sub_OldAvDist((cond-1)*2+loc)=length(curr_conf)/triads_per_cond; %accuracy rates
        for rate=1:num_ratings
            curr_sub_CountEachRatingAvDist((cond-1)*2*num_ratings+(loc-1)*num_ratings+rate)=sum(curr_conf==rate);
        end
        curr_conf=conf_btSTD(curr_trials,loc);
        curr_conf=curr_conf(curr_conf~=0); %remove if subject did not respond
        curr_sub_confSDAvDist((cond-1)*2+loc)=sum(curr_conf)/length(curr_conf);
    
    
    loc=2;
        curr_conf=conf_bt(curr_trials,loc)+conf_bt(curr_trials,(loc+1));
        curr_conf=curr_conf(curr_conf~=0); %take only those in which the participant chose this answer (loc=1 is the target)
        curr_sub_confAvDist((cond-1)*2+loc)=sum(curr_conf)/length(curr_conf); %average
        curr_sub_OldAvDist((cond-1)*2+loc)=length(curr_conf)/triads_per_cond; %accuracy rates
        for rate=1:num_ratings
            curr_sub_CountEachRatingAvDist((cond-1)*2*num_ratings+(loc-1)*num_ratings+rate)=sum(curr_conf==rate);
        end
        curr_conf=conf_btSTD(curr_trials,loc)+conf_btSTD(curr_trials,(loc+1));
        curr_conf=curr_conf(curr_conf~=0); %remove if subject did not respond
        curr_sub_confSDAvDist((cond-1)*2+loc)=sum(curr_conf)/length(curr_conf);
end


end
