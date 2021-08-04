function simPEL11_randomize_images(subj_num,subj_id)

%% initial definitions
rand('state',sum(100.*clock));%solve the randperm problem in matlab from being consistent
% [~, hostname]=system('hostname');
% if strcmp(hostname(1:6),'joanna')%the hostname command gives 1X7 char output, we only need the first 6.
%     project_dir='/Volumes/Oded/Bein/simPEL/simPEL11_onlyConsSameSimilarHalf';
% elseif strcmp(hostname(1:9),'electrode')%the hostname command gives 1X7 char output, we only need the first 6.
%     project_dir='/Volumes/Oded/Bein/simPEL/simPEL11_onlyConsSameSimilarHalf';
% else 
%     project_dir='/Volumes/davachilab/Bein/simPEL/simPEL11_onlyConsSameSimilarHalf';
% end

project_dir=pwd;
scripts_dir='scripts';
scritps_dir=fullfile(project_dir,scripts_dir);
subj_dir=fullfile(project_dir,'data',sprintf('%d%s',subj_num,subj_id));
if ~exist(subj_dir), mkdir(subj_dir); end
stim_dir='images';
stim_dir=fullfile(project_dir,stim_dir);
paired_dir='paired';

%cond: experimental conditions: violation/no-violation
%categ: size category of the stimuli: small, big
pairs_per_cond=36;
num_conds=2;%1-violation,2-no violation
categ_names={'small','big'};
num_categ=length(categ_names);

foils_per_categ=[pairs_per_cond pairs_per_cond]./2;%number of foils during item test in each category: small, big - since only half of the items will be presented, number of foils should correspond to that.
num_Ssmall=pairs_per_cond/num_categ; %number of single items small in each v/nv condition - make sure that pairs_per_cond/num_categ dvide to integers
num_Sbig=pairs_per_cond/num_categ; %number of single items big in each v/nv condition - make sure that pairs_per_cond/num_categ dvide to integers
single_per_categ=[num_Ssmall*num_conds num_Sbig*num_conds]; %total single items
items_per_cond_per_categ=pairs_per_cond*num_conds+foils_per_categ;%big/small items that wiil be sllocated to either foils/initial pairs in the study
reps_init_learning=12;

%get all the images - this is unnecessary, but good for debugging:
% loop through categories
for c=1:num_categ
    
    % get folder path for image dir
    tempStimDir=fullfile(stim_dir,categ_names{c});
    
    % get filenames
    dir_list{c} = dir(tempStimDir); %#ok<*AGROW>
    
    % skip stuff that isn't an image
    dir_list{c} = dir_list{c}(3:end);               % skip . & ..
    if (strcmp(dir_list{c}(1).name,'.DS_Store')==1) % also skip .DS_Store
        dir_list{c} = dir_list{c}(2:end);
    end
end

%get the paired images:
% loop through categories
for c=1:num_categ
    
    % get folder path for image dir
    tempStimDir=fullfile(stim_dir,paired_dir,categ_names{c});
    
    % get filenames
    paired_dir_list{c} = dir(tempStimDir); %#ok<*AGROW>
    
    % skip stuff that isn't an image
    paired_dir_list{c} = paired_dir_list{c}(3:end);               % skip . & ..
    if (strcmp(paired_dir_list{c}(1).name,'.DS_Store')==1) % also skip .DS_Store
        paired_dir_list{c} = paired_dir_list{c}(2:end);
    end
end
%%%%%prepare a file coupling the images with their numbers%%%%
%% setting the randomization of objects' sizes

%the response wanted - -1 means smaller 1 means bigger
%there will be big/small/medium objects in the experiment

%there are 36 orders in total: first 18 will be big-small, last 18 will be
%small-big
size_mat=[repmat([2,1],pairs_per_cond/2,1);repmat([1,2],pairs_per_cond/2,1)];
response_mat=[repmat([1,-1],pairs_per_cond/2,1);repmat([-1,1],pairs_per_cond/2,1)];
      
%now, for each condition, choose x items of each size, and assign them to
%pairs randomly.

%stim_alocation:
stim_rand=randperm(items_per_cond_per_categ(1))';
stim_rand=[stim_rand randperm(items_per_cond_per_categ(2))'];
    
%stim rand has two columns: small(10000),big(20000) each row is a pair
stim_rand=stim_rand+repmat([10000 20000],items_per_cond_per_categ(1),1);

%assign pairs to conditions: violation:1000, no_violation=2000
stim_v=stim_rand(1:pairs_per_cond,:)+1000;
stim_nv=stim_rand((pairs_per_cond+1):2*pairs_per_cond,:)+2000;

%in each condition, set the items in the correct order based on the size
%matrix
for row=1:size(size_mat,1)
    stim_v(row,:)=stim_v(row,size_mat(row,:));
    stim_nv(row,:)=stim_nv(row,size_mat(row,:));
end
%insert the wanted response into the sequence (1:bigger, -1: smaller)
stim_v=stim_v.*response_mat;
stim_nv=stim_nv.*response_mat;

%foils for test
stim_foils=[];
for c=1:num_categ
    stim_foils=[stim_foils; stim_rand((num_conds*pairs_per_cond+1):items_per_cond_per_categ(c),c)+9000];
end
stim_foils=stim_foils(randperm(length(stim_foils)));

%% grab items for all single items:2 conditions,in each there are Ssmall small items,
%Smed medium and Sbig big.

%first, choose which of the stimuli pairs will be in the study, and which
%will be lure, in each category.

%small:
single_stim_rand_temp=1:pairs_per_cond*2';
single_stim_rand_temp=reshape(single_stim_rand_temp,2,pairs_per_cond)';
%now randomize between the columns:
for i=1:pairs_per_cond
    single_stim_rand_temp(i,:)=single_stim_rand_temp(i,randperm(2));
end
single_stim_rand=single_stim_rand_temp(randperm(pairs_per_cond),:);

%big:
single_stim_rand_temp=1:pairs_per_cond*2';
single_stim_rand_temp=reshape(single_stim_rand_temp,2,pairs_per_cond)';
%now randomize between the columns:
for i=1:pairs_per_cond
    single_stim_rand_temp(i,:)=single_stim_rand_temp(i,randperm(2));
end
single_stim_rand=[single_stim_rand single_stim_rand_temp(randperm(pairs_per_cond),:)];

single_stim_rand=single_stim_rand+repmat([10000 10000 20000 20000],pairs_per_cond,1);

stim_single_v=[single_stim_rand(1:pairs_per_cond/2,1);single_stim_rand(1:pairs_per_cond/2,3)]+5000;%single v items targets
lures_single_v=[single_stim_rand(1:pairs_per_cond/2,2);single_stim_rand(1:pairs_per_cond/2,4)]+7000;%single v items lures only at test
stim_single_nv=[single_stim_rand((pairs_per_cond/2+1):pairs_per_cond,1);single_stim_rand((pairs_per_cond/2+1):pairs_per_cond,3)]+6000;%single nv items targets
lures_single_nv=[single_stim_rand((pairs_per_cond/2+1):pairs_per_cond,2);single_stim_rand((pairs_per_cond/2+1):pairs_per_cond,4)]+8000;%single nv items lures only at test
stim_single_v=stim_single_v.*response_mat(:,2);%insert the wanted response into the sequence (1:bigger, -1: smaller)
stim_single_nv=stim_single_nv.*response_mat(:,2);%insert the wanted response into the sequence (1:bigger, -1: smaller)

%evantually, each stimuli has a number:
%The ten thousands digit indicates the size: small(10000),big(20000)
%The thuosands digit indicates the condition: violation (to be):1000, no-violation (to be):2000, violation (in the violation block)=3000, violation (in the violation block)=4000,...
%                                             single violation=5000  single no-violation=6000 foil=9000.
%Then the stimul: 1-117 in each size.

%save it - need the same randomization for all parts
save([subj_dir '/items_randomization'],'stim_v','stim_nv','stim_single_v','stim_single_nv','lures_single_v','lures_single_nv','stim_foils');

%% construct trials sequences for initial learning sessions

%construct pair sequence for initial learning, day 1 (v,nv):
pairs_temp=[ones(1,pairs_per_cond)*1000 ones(1,pairs_per_cond)*2000];
pairs_temp=pairs_temp+[1:pairs_per_cond 1:pairs_per_cond];

pairs_seq_day1=zeros(reps_init_learning,length(pairs_temp));

counter=0;
good_seq=1;
while good_seq
    good_seq=0;
    counter=counter+1;
    for r=1:reps_init_learning
        check=1;
        while check %a loop to check 1-back or 2-back
            pairs_seq_day1(r,:)=pairs_temp(randperm(length(pairs_temp)));
            check=0;
            
            %make sure that there's no 1-back or 2-back identical pairs between
            %repetitions
            if r>=2 %strat checking only from second rep
                if pairs_seq_day1(r,1)==pairs_seq_day1(r-1,length(pairs_temp)) || pairs_seq_day1(r,1)==pairs_seq_day1(r-1,length(pairs_temp)-1)
                    check=1; %keep on shuffling
                elseif pairs_seq_day1(r,2)==pairs_seq_day1(r-1,length(pairs_temp)) %there is no risk that pairs_seq_day1(r,2)==pairs_seq_day1(r,1) so no need to check
                    check=1; %keep on shuffling
                end
            end
        end
    end
    
    if good_seq~=1 %the sequence is good so far, continue to allocate the stimuli
        %         counter
        %construct the stimuli sequence:
        stim_seq_day1=zeros(reps_init_learning,length(pairs_temp)*2);
        stim_loc_seq_day1=zeros(size(stim_seq_day1));
        stim_in_pair_seq_day1=zeros(size(stim_seq_day1));
        
        for r=1:reps_init_learning
            curr_loc=1;
            for t=1:length(pairs_temp)
                stim_loc_seq_day1(r,curr_loc:curr_loc+1)=[1 2];
                curr_pair_num=mod(pairs_seq_day1(r,t),1000);
                stim_in_pair_seq_day1(r,curr_loc:curr_loc+1)=curr_pair_num;
                if pairs_seq_day1(r,t)<2000 %this is a violation condition
                    stim_seq_day1(r,curr_loc:curr_loc+1)=stim_v(curr_pair_num,:);
                    curr_loc=curr_loc+2;
                elseif (pairs_seq_day1(r,t)>2000) %this is a no violation condition
                    stim_seq_day1(r,curr_loc:curr_loc+1)=stim_nv(curr_pair_num,:);
                    curr_loc=curr_loc+2;
                end
            end
        end
    end
    
    stim_in_pair_seq_day1(stim_in_pair_seq_day1==0)=pairs_per_cond;%correct the modulus
    
    if counter==100000
        break
    end
        
end %ends the while loop that create the sequence for the initial learning day 1

%% generate jitter initial learning - that will influence the randmoization of the reminder and the err block, so do it before.

%made it not jittered - all of that is irrelevant., here only for future
%ref - if we want to bring it back

%this will be out of the scanner. the reason that it is jittered is so that
%it could be jittered in the scanner as well. so, can keep the jittering
%identical across all repetitions of a sequence. The idea is that I could,
%in the scanner, have a regressor per condition:
%B items - violation, V items, C items - no violation, V items, violation.
%since this only occur once, jittering between these items will ensure that
%these univariate responses are separated. keeping the same jittering would
%mean that the B and C regressors in the no-violation conditions, if taking both repetitions, are correlated.
%but, since I'll always be interested in the reminder/violation condition -
%these will be modelled as 2 different regressors. within each repetition,
%it is jitterd.
%pos 1: btw 1st and 2nd items
%pos 2: btw 2nd and 3rd items
%pos 3: after 3rd item

%eventually, pos 2 will be pre-violation, and pos 3 will be post violation
%for the no-violation condition, it'll be pos 3 pre-no vioaltion, and, will
%need to match another column for the post-no-violation.
%so, pos 2 and pos 3 columns should be identical.
%the post no-violation column should be matched, such that for each gap in
%pos 3, find the corresponding in pos 2, and put the pos 3 that corresponds
%to it.

%%%%%later - will need to control the randmoization to have identical
%%%%%lengths of runs. that's important for the fMRI, not included here - so
%%%%%runs can be few seconds shorter/longer, doesn't seem critical to me (7.4.17).

big_j=perms(1:4);
% %isi=[1 1 3 5];

% 
% big_j=big_j(:,2:4);%need only 3 columns
% %there are 24 randomizations, and 32 pairs per condition, so
% %choose another 8: 2 of each combination of 3 numbers, I just chose the
% %ones that in the 2 position in big_j (to be pos 1), has 1,2,3,4
% %while doing that, organize in a manner that is easier to parse (for randmoization later):
% big_j=big_j([1:6 5:6 7:12 9:10 13:18 13:14 19:24 23:24],:);
% %now there are 8 in each combination of lengths:

%that is needed for the randomization of the pairs later on:
rand_pairs_isi=randperm(pairs_per_cond);
%two rows, dividing the pairs that will appear in the two sections of the reminder by similar isi combo, to equate the length
rand_pairs_rem=[rand_pairs_isi(1:2:end) rand_pairs_isi(2:2:end)]; 
% 
% isi_post_nv=zeros(size(big_j,1),1);
% %find the corresponding pos3 items to the pos2 items, and put them after
% %the pos3 items with the same gap - the idea is that across conditions, pre
% %and post v and nv gaps are equal - so that for every pre-post v pair,
% %threre is a corresponding pre-post nv pair
% %that's the randmoization method - but intractable for where is the
% %equivalent one in the isi_post_nv, so below I did it in a more tractable
% %way:
% for i=1:4
%     pos3_after2=big_j(big_j(:,2)==i,3); %find all the isi after 3 that correspond to the pos2 currently considered.
%     pos3_after2=pos3_after2(randperm(length(pos3_after2))); %shuffle them
%     isi_post_nv(big_j(:,3)==i)=pos3_after2; %put them after the corresponding pos3 items
% end
% %isi_post_3rd_item=repmat(isi,1,pairs_per_cond/length(isi))'; %note that number of total pairs must be divisable by 4 - otherwise, assign the values differently.
% % isi_post_nv=[2:2:pairs_per_cond;1:2:pairs_per_cond];
% % isi_post_nv=reshape(isi_post_nv,pairs_per_cond,1);
% % isi_post_nv=big_j(isi_post_nv,3);
% % big_j_post_nv=[];
% % for i=1:4
% %     big_j_post_nv=[big_j_post_nv;big_j(isi_post_nv==i,:) isi_post_nv(isi_post_nv==i)];
% % end
% % %big_j_post_nv=[big_j_post_nv isi_post_nv];
% % big_j_post_nv=[big_j isi_post_nv];
% 
% isi_day1=zeros(size(stim_seq_day1));
% for r=1:reps_init_learning
%     for t=1:pairs_per_cond
%         %should be two pairs in each rep
%         isi_day1(r,stim_in_pair_seq_day1(r,:)==rand_pairs_isi(t))=isi(repmat(big_j(t,:),1,num_conds)); %row t in big_j will go to pair (rand_pairs_isi(t))
%     end
% end
% 


%parameters of the simulation
run_delay = 2;
run_decay = 2;
stim_dur = 1.5;
isi=0.5;
isi_day1=ones(size(stim_seq_day1))*isi;
trial_lengths=stim_dur + isi_day1;
day1_init_onsets=run_delay + cumsum(trial_lengths,2) - (stim_dur+isi_day1);


%% construct trials sequences for the reminder block
runs_rem=2; %divided to two runs
%construct pair squence for reminder block on day 2 (1000: v reminder, 2000:nv reminder):
pairs_temp=[ones(1,pairs_per_cond)*1000 ones(1,pairs_per_cond)*2000];

pairs_temp=pairs_temp+[rand_pairs_rem rand_pairs_rem]; %take the pairs based on the isi randomization
seq_length=length(pairs_temp);
%make sure that in each half there is the same number of v and nv pairs:
pairs_seq_reminder=[pairs_temp([1:seq_length/4 seq_length/2+1:seq_length*3/4]);...
                     pairs_temp([seq_length/4+1:seq_length/2 seq_length*3/4+1:seq_length])];

for r=1:runs_rem
    pairs_seq_reminder(r,:)=pairs_seq_reminder(r,randperm(seq_length/runs_rem));
end

%construct the stimuli sequence:
stim_seq_reminder=zeros(runs_rem,seq_length/runs_rem*2);
stim_loc_seq_reminder=zeros(size(stim_seq_reminder));
stim_in_pair_seq_reminder=zeros(size(stim_seq_reminder));

for r=1:runs_rem
    curr_loc=1;
    for t=1:seq_length/2
        curr_pair_num=mod(pairs_seq_reminder(r,t),1000);
        curr_cond=floor(pairs_seq_reminder(r,t)/1000);
        stim_loc_seq_reminder(r,curr_loc:curr_loc+1)=[1 2];
        stim_in_pair_seq_reminder(r,curr_loc:curr_loc+1)=curr_pair_num;
        switch curr_cond
            case 1 %this is a violation condition
                stim_seq_reminder(r,curr_loc:curr_loc+1)=stim_v(curr_pair_num,:);
                curr_loc=curr_loc+2;
            case 2 %this is a no-violation condition
                stim_seq_reminder(r,curr_loc:curr_loc+1)=stim_nv(curr_pair_num,:);
                curr_loc=curr_loc+2;
        end
    end
end
stim_in_pair_seq_reminder(stim_in_pair_seq_reminder==0)=pairs_per_cond;%correct the reminder

%% generate jitter reminder sessions - don't check correlations (but has set up here, if needed)
% run_length=364;
% run_length_tp=run_length/tr;
% each repetition is now on a different scan
% check=1;
% while check
%     check=0;


%parameters of the simulation
run_delay = 2;
run_decay = 2;
tr=2;

% % stick function convolved with FSL's double gamma HRF
hrf =[0.0036 0.0960 0.1469 0.1163 0.0679 0.0329 0.0137 0.0047 0.0008 -0.0007 -0.0013 -0.0017 -0.0017 -0.0017 -0.0017 -0.0017 -0.0017];
%
%% put in the jitter:

%no jitter:
isi_reminder=ones(size(stim_seq_reminder))*isi;

trial_lengths=stim_dur + isi_reminder;
reminder_onsets=run_delay + cumsum(trial_lengths,2) - (stim_dur+isi_reminder);

% isi_reminder=zeros(size(stim_seq_reminder'));
% for t=1:pairs_per_cond
%     isi_reminder(stim_in_pair_seq_reminder'==rand_pairs_isi(t))=isi(repmat(big_j(t,:),1,num_conds)); %find the equivalent pairs in all conditions and put the isi in all of them
%      %isi_day1(r,stim_in_pair_seq_day1(r,:)==rand_pairs_isi(t))=isi(repmat(big_j(t,:),1,num_conds)); %row t in big_j will go to pair (rand_pairs_isi(t))
% end

% %check run length
% rem_run_lengths=zeros(runs_rem,1);
% for i=1:runs_rem
%     rem_run_lengths(i)=reminder_onsets(i,end)+isi_reminder(i,end);
% end
% 
% if any(rem_run_lengths~=rem_run_lengths(1))
%     fprintf('script is bad, unequal reminder run lengths\n');
% end

% 
%     %now check the correlations, for each run:
%     for r=1:runs_rem
%             
%             trial_onsets_tp=reminder_onsets./tr;
%             design=zeros(run_length_tp,6);
%             design_conv=[];
%             %now check for each position: A,B,C, v/nv:
%             for cond=1:2
%                 for loc=1:3
%                     col=(cond-1)*3+loc;
%                     curr_trials=floor(mod(abs(stim_seq_reminder(r,:)),10000)/1000)==cond;
%                     curr_trials(stim_loc_seq_reminder(r,:)~=loc)=0; %zero all those that are not in the current location
%                     design(trial_onsets_tp(curr_trials),col)=1;
%                     design_conv(:,col)=conv(design(:,col),hrf);
%                 end
%             end
%             design_conv=design_conv(1:run_length_tp,:);
%             design_corr=abs(corr(design_conv));
%             design_corr_vec=triu(design_corr)~=0; %CHECK THIS
%             if any(any(design_corr_vec>min_corr))
%                 jitter_check=1;
%             end
% 
%     end
%end

%% construct trials sequences for reminder/error-induction block

runs_violation=4; %divide to 4 scans, each will have 8 pairs in each condition, reminder and v/nv
pairs_per_run=pairs_per_cond/runs_violation;
%construct pair squence for reminder/violation runs, day 1 (1000: v reminder, 2000:nv reminder).
%randomize the reminder, the violation will be set by that (see below). that creates a randomization in which in each run there are 8 pairs in
%each condition. Each will appear twice: reminder and v/nv :

%after randomizing - making sure that the run lengths are equal - no need
%here because we have fixed isi, but if needed for the furture.
check_length=1;
attempts_length=0;
%vnv_corr_length=426.5;%in secs - that's what it should be. I ran a few simulations and checked that.
while check_length
    check_length=0;
    attempts_length=attempts_length+1;
    pairs_temp=[ones(pairs_per_cond,1)*1000 ones(pairs_per_cond,1)*2000];
    rand_v=randperm(pairs_per_cond);
    rand_nv=randperm(pairs_per_cond);
    pairs_temp=pairs_temp+[rand_v' rand_nv'];
    pairs_temp=reshape(pairs_temp,runs_violation,pairs_per_cond*num_conds/runs_violation);
    pairs_v_rem=pairs_temp(:,1:pairs_per_run);
    pairs_nv_rem=pairs_temp(:,pairs_per_run+1:pairs_per_run*2);
    %
    
    pairs_seq_vnv=zeros(runs_violation,pairs_per_run*num_conds*2);
    attempts=1;
    rand_conds=[ones(1,5) ones(1,5)*2 ones(1,6)*3 ones(1,6)*4 ones(1,3)*3 ones(1,3)*4]; %9 pairs per run in each condition: the begining is only reminder, so set elsewhere, then: 5 more reminder trials in each condition, 9 violation trials, the last 6 will always be v/nv trials - because there is a minimun of 6 pairs btw reminber and v/nv
    for r=1:runs_violation
        %fprintf(sprintf('prep run %d vnv block\n',r));
        check=1;
        while check %a loop to check 1-back or 2-back
            check=0;
            temp_seq_rem=[pairs_v_rem(r,1:4) pairs_nv_rem(r,1:4) pairs_v_rem(r,5:9) pairs_nv_rem(r,5:9)];
            temp_seq_rem=temp_seq_rem([randperm(8) randperm(10)+8]); %start with 4 of each condition, then randomize
            %to make sure that the violation always appear after the reminder,
            %shuffle within each 4, then within each 5 (total of 18 pairs) and place the quartet one after the other to construct the sequence
            temp_seq_viol=[temp_seq_rem(randperm(4))+2000 temp_seq_rem(randperm(4)+4)+2000];
            temp_seq_viol=[temp_seq_viol temp_seq_rem(randperm(5)+8)+2000 temp_seq_rem(randperm(5)+13)+2000];

            temp_seq=zeros(1,pairs_per_run*4);
            temp_seq(1:8)=temp_seq_rem(1:8);
            temp_seq_rem=temp_seq_rem(9:end);
            rand_conds=rand_conds([randperm(length(rand_conds)-6) (randperm(6)+length(rand_conds)-6)]); %randomize separately all the listminus the last 6, and the last six items.
            for i=1:2
                temp_seq(find(rand_conds==i)+8)=temp_seq_rem(floor(temp_seq_rem/1000)==i);
            end
            for i=3:4
                temp_seq(find(rand_conds==i)+8)=temp_seq_viol(floor(temp_seq_viol/1000)==i);
            end
            
            %now check that the sequence is good:
            
            attempts=attempts+1;
            %make sure that there are no more than 3 consequtive error pairs
            div_pairs=double(floor(temp_seq/1000)==3); %mark the violation - 3
            find_conseq_er=conv(div_pairs,[0.25 0.25 0.25 0.25]); %average across 4 pairs
            if any(find_conseq_er(4:length(div_pairs))==1) %if the average across 4 pairs is 1, that means there are 4 consecutive error pairs - no good.
                check=1;
            else %continue to check that all violations are at least 7 pairs (6 pairs gap) after their reminder
                %fprintf('checking the gaps btw violations and reminder run %d\n',r);
                check_viol=temp_seq(temp_seq>3000);
                for i=1:length(check_viol)
                    loc_viol=find(temp_seq==check_viol(i));
                    loc_rem=find(temp_seq==(check_viol(i)-2000));
                    if (loc_viol-loc_rem<7) %|| (loc_viol-loc_rem>15)
                        check=1;
                        break
                    end
                end
            end
            
            if attempts>1000000
                fprintf('didn''t find a good violation run sequence\n');
                check=1;
            end
            
        end
        if ~check %found a good sequence:
            pairs_seq_vnv(r,:)=temp_seq;
        end
        
    end
    
    %% another way of checking the distance between pairs
    %         if ~check %sequence is good so far
    %             %check that there is a minimum distance of 6 pairs between the
    %             %reminder and the violation/no-violation
    %             for i=9:length(temp_seq) %only check from the 9th item - I set up the first 8 to be reminder
    %                 curr_check=temp_seq(i-5:i);
    %                 diff=repmat(curr_check,length(curr_check),1)-repmat(curr_check',1,length(curr_check));%compute all differences
    %                 if any(diff==2000) %anything is exactly 2000 - that means that there is a close pair
    %                     check=1;
    %                     break
    %                 end
    %             end
    %         end
    %% construct the stimuli sequence,and insert the single items:
    
    %(violation: 32) - will be inseated instead of third items.
    %(no violation: 32) - will be insearted after the nv pair
    stim_seq_vnv=zeros(runs_violation,size(pairs_seq_vnv,2)*2+pairs_per_cond/runs_violation); %in all conditions there are 2 items per pair, and add pairs_per_cond/runs_violation more for the nv items that are added to nv pairs
    stim_loc_seq_vnv=zeros(size(stim_seq_vnv)); 
    stim_in_pair_seq_vnv=zeros(size(stim_seq_vnv));
    
    %%%% for the nv single items - they are not insearted instead of the C item,
    % but after it. So to maintain the response identical to the violation condition, they should be shuffled, to
    % match, so that big items will appear after small items, and small
    % items after big items
    rand_nv_single=zeros(pairs_per_cond,1);
    rand_nv_single(1:num_Ssmall)=randperm(num_Ssmall)+num_Sbig;
    rand_nv_single((num_Ssmall+1):end)=randperm(num_Ssmall);
    
    for r=1:runs_violation
        curr_loc=1;
        for t=1:size(pairs_seq_vnv,2)
            curr_pair_num=mod(pairs_seq_vnv(r,t),1000);
            curr_cond=floor(pairs_seq_vnv(r,t)/1000);
            stim_loc_seq_vnv(r,curr_loc:curr_loc+1)=[1 2];
            stim_in_pair_seq_vnv(r,curr_loc:curr_loc+1)=curr_pair_num;
            switch curr_cond
                case 1 %this is a violation condition - reminder - just put the sequence
                    stim_seq_vnv(r,curr_loc:curr_loc+1)=stim_v(curr_pair_num,1:2); %first two items remain the same
                    curr_loc=curr_loc+2;
                    
                case 2 %this is a no violation condition - reminder - just put the sequence
                    stim_seq_vnv(r,curr_loc:curr_loc+1)=stim_nv(curr_pair_num,1:2);
                    curr_loc=curr_loc+2;
                case 3 %this is a violation condition - violation
                    stim_seq_vnv(r,curr_loc)=stim_v(curr_pair_num,1); %first item remains the same
                    %mark them as violation round - by making the thousands digits 3
                    for cc_loc=curr_loc
                        if stim_seq_vnv(r,cc_loc) < 0
                            stim_seq_vnv(r,cc_loc)=stim_seq_vnv(r,cc_loc)-2000;
                        else
                            stim_seq_vnv(r,cc_loc)=stim_seq_vnv(r,cc_loc)+2000;
                        end
                    end
                    curr_loc=curr_loc+1;
                    %inseart the proper single item:
                    stim_seq_vnv(r,curr_loc)=stim_single_v(curr_pair_num); %single item x apears after pair x
                    curr_loc=curr_loc+1;
                case 4 %this is a no violation condition - no violation - put the sequence, and add the item
                    stim_seq_vnv(r,curr_loc:curr_loc+1)=stim_nv(curr_pair_num,:); %all items remain the same'
                    %mark them as no violation round - by making the thousands digits 4
                    for cc_loc=curr_loc:curr_loc+1
                        if stim_seq_vnv(r,cc_loc) < 0
                            stim_seq_vnv(r,cc_loc)=stim_seq_vnv(r,cc_loc)-2000;
                        else
                            stim_seq_vnv(r,cc_loc)=stim_seq_vnv(r,cc_loc)+2000;
                        end
                    end
                    curr_loc=curr_loc+2;
                    %inseart the proper single item:
                    stim_seq_vnv(r,curr_loc)=stim_single_nv(rand_nv_single(curr_pair_num)); %single item rand_nv_single(x) appears after pair x
                    stim_in_pair_seq_vnv(r,curr_loc)=curr_pair_num;
                    stim_loc_seq_vnv(r,curr_loc)=4;
                    curr_loc=curr_loc+1;
            end
        end
    end
    
    %check the sequence:
    if length(unique(abs(stim_seq_vnv))) ~= (pairs_per_cond*num_conds*2+pairs_per_cond*2+pairs_per_cond*3)%reminder + violation pairs + nv pairs
        fprintf('not showing all items, and each one only once, in stim_seq_vnv\n');
    end
    
    stim_seq_vnvCheck=stim_seq_vnv'; %transpose just so it'll be easier to check
    for i=1:pairs_per_cond
        aItem=stim_v(i,1);
        if aItem < 0
            aItem=aItem-2000;
        else
            aItem=aItem+2000;
        end
        aItem_viol=find(stim_seq_vnvCheck==aItem);
        vItem=stim_seq_vnvCheck(aItem_viol+1);
        vItem_inStims=find(stim_single_v==vItem);
        %fprintf(sprintf('check violation pair %d\n',i));
        if vItem_inStims ~= i
            fprintf('a violation item is not after its matching A item\n');
        end
    end
    
    for i=1:pairs_per_cond
        bItem=stim_nv(i,2);
        if bItem < 0
            bItem=bItem-2000;
        else
            bItem=bItem+2000;
        end
        bItem_viol=find(stim_seq_vnvCheck==bItem);
        vItem=stim_seq_vnvCheck(bItem_viol+1);
        vItem_inStims=find(stim_single_nv==vItem);
        %fprintf(sprintf('check no violation pair %d\n',i));
        if vItem_inStims ~= rand_nv_single(i)
            fprintf('a no violation item is not after its matching B item\n');
        end
    end
    
    stim_in_pair_seq_vnv(stim_in_pair_seq_vnv==0)=pairs_per_cond;%correct the reminder
    
    %% generate jitter error sessions - and check the gap
    isi_vnv_block=ones(size(stim_seq_vnv))*isi;
    trial_lengths=stim_dur + isi_vnv_block;
    vnv_onsets=run_delay + cumsum(trial_lengths,2) - (stim_dur+isi_vnv_block);

    
%     %here its redundent since it's all the same isi, but I kept it for the
%     %future
%     %parameters of the simulation
%     run_delay = 12;
%     run_decay = 12;
%     tr=2;
%     
%     % % stick function convolved with FSL's double gamma HRF
%     hrf =[0.0036 0.0960 0.1469 0.1163 0.0679 0.0329 0.0137 0.0047 0.0008 -0.0007 -0.0013 -0.0017 -0.0017 -0.0017 -0.0017 -0.0017 -0.0017];
%     %
%     %each repetition is now on a different scan
%     %DON'T CHECK CORRELATIONS IN THIS VERSION
%     isi_vnv_block=zeros(size(stim_seq_vnv'));
%     stim_loc_seq_vnv=stim_loc_seq_vnv';
%     for t=1:pairs_per_cond
%         pairs_loc=find(stim_in_pair_seq_vnv'==rand_pairs_isi(t));
%         pair_pos=stim_loc_seq_vnv(pairs_loc);
%         %isi_reminder(stim_in_pair_seq_reminder'==rand_pairs_isi(t))=isi(repmat(big_j(t,:),1,num_conds));
%         isi_vnv_block(pairs_loc(pair_pos~=4))=isi(repmat(big_j(t,:),1,num_conds*2)); %twice for reminder and v/nv
%         isi_vnv_block(pairs_loc(pair_pos==4))=isi(isi_post_nv(t));
%     end
%     isi_vnv_block=isi_vnv_block';
%     stim_loc_seq_vnv=stim_loc_seq_vnv';
%     trial_lengths=stim_dur + isi_vnv_block;
%     vnv_onsets=run_delay + cumsum(trial_lengths,2) - (stim_dur+isi_vnv_block);
%     
%     %now check the lengths
%     vnv_run_lengths=zeros(runs_violation,1);
%     for i=1:runs_violation
%         %fprintf('computing violation seq lengths\n');
%         vnv_run_lengths(i)=vnv_onsets(i,end)+stim_dur+isi_vnv_block(i,end);
%     end
%     %fprintf('checking violation seq lengths\n');
%     if any(vnv_run_lengths~=vnv_corr_length)
%         %fprintf('violation seq lengths are not good\n');
%        check_length=1; 
%     end
%     
%     if ~check_length

%check the gap between the violations
gap_v=[];
gap_nv=[];
cond_vnv_check=abs(stim_seq_vnvCheck);
cond_vnv_check=floor(rem(cond_vnv_check,10000)/1000);
vnv_onsetsCheck=vnv_onsets';
for t=1:pairs_per_cond
    pair_loc=find(((stim_in_pair_seq_vnv'==t)+(cond_vnv_check==1))==2,1); %find the A item - rem
    pair1_t=vnv_onsetsCheck(pair_loc);
    pair_loc=find(((stim_in_pair_seq_vnv'==t)+(cond_vnv_check==3))==2,1); %find the A item - violation
    pair2_t=vnv_onsetsCheck(pair_loc);
    if (pair2_t-pair1_t)< 0
        fprintf('(pair2_t-pair1_t)< 0 \n');
    end
    gap_v=[gap_v;(pair2_t-pair1_t)];
end

for t=1:pairs_per_cond
    pair_loc=find(((stim_in_pair_seq_vnv'==t)+(cond_vnv_check==2))==2,1); %find the A item - rem
    pair1_t=vnv_onsetsCheck(pair_loc);
    pair_loc=find(((stim_in_pair_seq_vnv'==t)+(cond_vnv_check==4))==2,1); %find the A item - nv
    pair2_t=vnv_onsetsCheck(pair_loc);
    if (pair2_t-pair1_t)< 0
        fprintf('(pair2_t-pair1_t)< 0 \n');
    end
    gap_nv=[gap_nv;(pair2_t-pair1_t)];
end

if abs(mean(gap_v)-mean(gap_nv))>1
   % fprintf('mean(gap_v)>mean(gap_nv) \n');
    check_length=1;
end

if attempts_length>100000
    fprintf('run script again, didn''t find vnv sequences of equal length\n');
    check_length=0;
end

    %now check the correlations, for each run:
    %I checked A, B, C, (separated first and second time) and violation and no-violation
    %also: model all A together (first time), A together (second time) and B
    %and C (first time), and B and C (second time). - predictable vs.
    %unpredictable items
    % for r=1:runs_violation
    %     run_length=vnv_onsets(r,end)+10;
    %     run_length_tp=run_length/tr;
    %     col=1;
    %     trial_onsets_tp=vnv_onsets./tr;
    %     design=zeros(run_length_tp,4*4);
    %     design_conv=[];
    %     %now check for each position: A,B,C, v/nv:
    %     for cond=1:4
    %         for loc=1:3
    %             if  ~(cond==3 && loc==3) %skip the C items in the error condition
    %             curr_trials=floor(mod(abs(stim_seq_vnv(r,:)),10000)/1000)==cond;
    %             curr_trials(stim_loc_seq_vnv(r,:)~=loc)=0; %zero all those that are not in the current location
    %             design(trial_onsets_tp(r,curr_trials),col)=1;
    %             design_conv(:,col)=conv(design(:,col),hrf);
    %             col=col+1;
    %             end
    %         end
    %     end
    %     for cond=5:6
    %         curr_trials=floor(mod(abs(stim_seq_vnv(r,:)),10000)/1000)==cond;
    %         design(trial_onsets_tp(r,curr_trials),col)=1;
    %         design_conv(:,col)=conv(design(:,col),hrf);
    %         col=col+1;
    %     end
    %     %now combine conditions
    %     %all A together (first rep):
    %     curr_trials=floor(mod(abs(stim_seq_vnv(r,:)),10000)/1000)<3;
    %     curr_trials(stim_loc_seq_vnv(r,:)~=1)=0; %zero all those that are not in the first location
    %     design(trial_onsets_tp(r,curr_trials),col)=1;
    %     design_conv(:,col)=conv(design(:,col),hrf);
    %     col=col+1;
    %      %all B and C together (first rep):
    %     curr_trials=floor(mod(abs(stim_seq_vnv(r,:)),10000)/1000)<3;
    %     curr_trials(stim_loc_seq_vnv(r,:)==1)=0; %zero all those that are in the first location
    %     design(trial_onsets_tp(r,curr_trials),col)=1;
    %     design_conv(:,col)=conv(design(:,col),hrf);
    %     col=col+1;
    %      %all A together (second rep):
    %     curr_trials=floor(mod(abs(stim_seq_vnv(r,:)),10000)/1000)==3;
    %     curr_trials=curr_trials+floor(mod(abs(stim_seq_vnv(r,:)),10000)/1000)==4;
    %     curr_trials(stim_loc_seq_vnv(r,:)~=1)=0; %zero all those that are not in the first location
    %     design(trial_onsets_tp(r,curr_trials),col)=1;
    %     design_conv(:,col)=conv(design(:,col),hrf);
    %     col=col+1;
    %      %all B and C together (second rep):
    %     curr_trials=floor(mod(abs(stim_seq_vnv(r,:)),10000)/1000)==3;
    %     curr_trials=curr_trials+floor(mod(abs(stim_seq_vnv(r,:)),10000)/1000)==4;
    %     curr_trials(stim_loc_seq_vnv(r,:)==1)=0; %zero all those that are in the first location
    %     design(trial_onsets_tp(r,curr_trials),col)=1;
    %     design_conv(:,col)=conv(design(:,col),hrf);
    %
    %     design_conv=design_conv(1:run_length_tp,:);
    %     design_corr=abs(corr(design_conv));
    % %     design_corr_vec=triu(design_corr)~=0; %CHECK THIS
    % %     if any(any(design_corr_vec>min_corr))
    % %         jitter_check=1;
    % %     end
    %
    % end
end


vnv_gap=mean(gap_v)-mean(gap_nv);
fprintf(sprintf('gap_v-gap_nv is %d \n',vnv_gap));
 
        
%prepare the onsets array for the reminder and error sections
save([subj_dir '/trial_sequences.mat'],'stim_seq_day1','stim_seq_reminder','stim_seq_vnv',...
                                       'isi_day1','isi_reminder','isi_vnv_block',...
                                       'stim_loc_seq_day1','stim_in_pair_seq_day1','day1_init_onsets',...
                                       'stim_loc_seq_reminder','stim_in_pair_seq_reminder','reminder_onsets',...
                                       'stim_loc_seq_vnv','stim_in_pair_seq_vnv','vnv_onsets','vnv_gap');

%% prepare trial sequence for the recognition test 

%the same/similar/lure randomization: in each of the conditions, half of
%the items would appear first as same, and half first as lures. because I
%wanted to have that, I didn't constrain the order so that it'll not be two
%constrained, or that the same items that apear as same (e.g.), would apear
%shortly after as lures

%create a pseudo-random sequence of trials, in which I control that each
%quarter of the test will have the same amount of trials from each
%condition (not to have a situation of all the pairs from one conditions
%being at the end of the test - importent when there's so little pairs

nRuns=4; %just to have the cond allocation correct - divide to half
all_conds=[ones(1,(pairs_per_cond*2/nRuns-1))*9 ones(1,pairs_per_cond/nRuns)*5 ones(1,pairs_per_cond/nRuns)*6 ones(1,pairs_per_cond/nRuns)*7 ones(1,pairs_per_cond/nRuns)*8]; %new, v-same and nv-same, v-lure, nv-lure items
%cond_seq=zeros(nRuns,length(all_conds));

nRuns=2;%actual number of runs
cond_seq=[];
for r=1:nRuns
    check=1;
    while check==1
        check=0;%assume that the randomization is good, to be changed if checks fail
        
        %for each quarter, leave the first one to be new (to equate for the
        %pre-trial durations:
        cond_seq_temp=[9 all_conds(randperm(length(all_conds)))];

        %check that there are no more than 4 identical responses in a row.
        for t=5:length(cond_seq_temp) %
            if all((cond_seq_temp(t-4:t)==5) + (cond_seq_temp(t-4:t)==6)) %"same" response
                check=1;
                break;
            elseif all((cond_seq_temp(t-4:t)==7) + (cond_seq_temp(t-4:t)==8))%"similar" response
                check=1;
                break;
            elseif all(cond_seq_temp(t-4:t)==9) %"new" response
                check=1;
                break;
            end
        end
    end
    cond_seq=[cond_seq; cond_seq_temp];
end

%% put items in the slots
%in simPEL6 or 7 - i did the randomization considering the blocks in the vnv block in this way:  
%alocate the v and nv trials - first half - coordinate the randomization
%with the vnv block, first half in the vnv block will be tested first, in a
%random order, then, same thing for the second half.

%Here, I don't take blocks into account because that limits the
%randomization due to the same/similar, so look at the other script if want
%to bring back randmoization taking vnv block into account.
cond_seq=cond_seq'; %for placing the items and checking them
min_gap=10; %make sure there is a gap of at least 10 trials between an item and it's lure

check=1;
attempts=0;

while check==1
    attempts=attempts+1;
    check=0;%assume that the randomization is good, to be changed if checks fail
    
    %randomize:
    rand_v_same=randperm(pairs_per_cond)';
    rand_v_lure=rand_v_same((pairs_per_cond/2+1):end); %take the second half of the randomization
    rand_v_lure=rand_v_lure(randperm(pairs_per_cond/2));%randmoize the order
    rand_v_same=rand_v_same(1:pairs_per_cond/2); %cut that to half
    
    rand_nv_same=randperm(pairs_per_cond)';
    rand_nv_lure=rand_nv_same((pairs_per_cond/2+1):end); %take the second half of the randomization
    rand_nv_lure=rand_nv_lure(randperm(pairs_per_cond/2));%randmoize the order
    rand_nv_same=rand_nv_same(1:pairs_per_cond/2); %cut that to half
    
    %place all items
    temp_trial_seq=zeros(size(cond_seq));
    temp_trial_seq(cond_seq==9)=stim_foils;
    
    temp_trial_seq(cond_seq==5)=stim_single_v(rand_v_same);
    temp_trial_seq(cond_seq==6)=stim_single_nv(rand_nv_same);
    temp_trial_seq(cond_seq==7)=lures_single_v(rand_v_lure);
    temp_trial_seq(cond_seq==8)=lures_single_nv(rand_nv_lure);
    
    %check the gaps - no need becuase not presenting the same items.
    %see simPEL8 scripts if needed!
    
    if mod(attempts,10000)==0
        attempts
    end
    
    if attempts==500000
        check=0;
        fprintf('re-run the script, didn''t find a good sequence for the recognition test \n');
    end

end

%save it
rec_trial_seq=temp_trial_seq';
cond_seq=cond_seq';


% 
% %check the allocation
% half_seq=[find(cond_seq(1:(length(cond_seq)/2))==5) find(cond_seq(1:(length(cond_seq)/2))==6)];
% half_seq=rec_trial_seq(half_seq);
% half_vnv_seq=stim_seq_vnv(1:2,:);
% half_vnv_seq=[half_vnv_seq(cond_vnv_check(1:2,:)==5) ;half_vnv_seq(cond_vnv_check(1:2,:)==6)];
% if any(sort(half_seq')~=sort(half_vnv_seq))
%     fprintf('first half rec does not correspond to first half, vnv\n');
% end
% half_seq=[find(cond_seq(((length(cond_seq)/2)+1):end)==5) find(cond_seq(((length(cond_seq)/2)+1):end)==6)];
% half_seq=half_seq+length(cond_seq)/2;
% half_seq=rec_trial_seq(half_seq);
% half_vnv_seq=stim_seq_vnv(3:4,:);
% half_vnv_seq=[half_vnv_seq(cond_vnv_check(3:4,:)==5) ;half_vnv_seq(cond_vnv_check(3:4,:)==6)];
% if any(sort(half_seq')~=sort(half_vnv_seq))
%      fprintf('second half rec does not correspond to second half, vnv\n');
% end

%%
%used to be jittered, not any more, but cod is as if it is jittered... :)
isi=ones(1,9)*3;
isi_rec=zeros(size(rec_trial_seq));
for r=1:nRuns
    for t=5:8
        isi_rec(r,find(cond_seq(r,:)==t)-1)=isi(randperm(length(isi)));
    end
    t=9;
    curr_loc=find(cond_seq(r,:)==t);
    %remove the first one:
    curr_loc=curr_loc(2:end);
    isi_rec(r,curr_loc-1)=[isi(randperm(length(isi)-1)+1) isi(randperm(length(isi)))];
end

%parameters of the simulation - don't simulate now
run_delay = 2;
run_decay = 2;
stim_dur = 3;

trial_lengths=stim_dur + isi_rec;
rec_onsets=run_delay + cumsum(trial_lengths,2) - (stim_dur+isi_rec);

% trial_lengths=stim_dur + isi_day1;
save([subj_dir '/rec_trial_sequences.mat'],'rec_trial_seq','rand_v_same','rand_nv_same','rand_v_lure','rand_nv_lure','rec_onsets');

%% explicit memory test - PART 1: testing for A-X (v), B-C (nv)
%These are the distractors for the explicit memory test. since I don't test
%for A or B items in the item test, no worries about that
smallDist=[[3:(pairs_per_cond/2) 1:2]' [4:(pairs_per_cond/2) 1:3]'];
bigDist=smallDist+pairs_per_cond/2;
AllDist=[smallDist;bigDist];

All_stim(:,:,1)=stim_v;
All_stim(:,:,2)=stim_nv;
num_conds=size(All_stim,3);

ExpTest_cues=squeeze(All_stim(:,1,:)); %all A items
%but actually, for the nv, we want the B items
ExpTest_cues(:,2)=All_stim(:,2,2);

ExpTest_dist=zeros(pairs_per_cond,3,num_conds); %each row is a pair, column 1: v/nv target. column 2: dist. coumn 3: dist. Each z level is a condition
ExpTest_dist(:,1,1)=stim_single_v(1:pairs_per_cond);%grab all single,v
ExpTest_dist(:,1,2)=stim_single_nv(rand_nv_single);%grab all single, nv

%check that the target is always correct
for i=1:pairs_per_cond
    aItem=ExpTest_cues(i,1);
    if aItem < 0
        aItem=aItem-2000;
    else
        aItem=aItem+2000;
    end
    trgt=ExpTest_dist(i,1,1);
    aItem_viol=find(stim_seq_vnvCheck==aItem);
    vItem=stim_seq_vnvCheck(aItem_viol+1);
    %fprintf(sprintf('check explicit part1 violation pair  %d\n',i));
    if vItem ~= trgt 
       fprintf('incorrect target, violation part 1\n');
    end
end

for i=1:pairs_per_cond
    bItem=ExpTest_cues(i,2);
    if bItem < 0
        bItem=bItem-2000;
    else
        bItem=bItem+2000;
    end
    trgt=ExpTest_dist(i,1,2);
    bItem_viol=find(stim_seq_vnvCheck==bItem);
    vItem=stim_seq_vnvCheck(bItem_viol+1);
    %fprintf(sprintf('check explicit part1 no violation pair  %d\n',i));
    if vItem ~= trgt 
       fprintf('incorrect target, violation part 1\n');
    end
end

ExpTest_dist(:,2:3,1)=[stim_single_v(AllDist(:,1)) stim_single_v(AllDist(:,2))];%grab all single_v items
ExpTest_dist(:,2:3,2)=[stim_single_nv(AllDist(rand_nv_single,1)) stim_single_nv(AllDist(rand_nv_single,2))];%grab all single_nv

%build the locations matrix:
loc=perms([1,2,3]);
all_locs=repmat(loc,3,1);%to have 18 - for big and small.
loc_allConds=zeros(pairs_per_cond,3,num_conds);
for cond=1:num_conds %each condition
    loc_allConds(:,:,cond)=[all_locs(randperm(size(all_locs,1)),:);all_locs(randperm(size(all_locs,1)),:)];
end

%build the stimuli in according to their locations
%so loc_allConds - if 3 is on the first colum - it means that item number 3
%is on the first location. so to look for the targets - it's where "1"s are
for cond=1:num_conds
    for pair=1:pairs_per_cond
    ExpTest_dist_inLoc(pair,:,cond)=ExpTest_dist(pair,loc_allConds(pair,:,cond),cond);
    end
end

%construct the order of the pairs - make sure that there's a gap between identical items;
%first three columns are the distractors, the cue is the forth column.
check=1;
attempts=1;
while check
    check=0;
    attempts=attempts+1;
    ExpTest_all_trials_InLoc=[];
    loc_allCondsInOrder=[];
    for cond=1:num_conds
        ExpTest_all_trials_InLoc=[ExpTest_all_trials_InLoc; [ExpTest_dist_inLoc(:,:,cond) ExpTest_cues(:,cond)]];
        loc_allCondsInOrder=[loc_allCondsInOrder;loc_allConds(:,:,cond)];
    end
    %shuffle:
    rand_parameter=randperm(length(ExpTest_all_trials_InLoc))';
    ExpTest_all_trials_InLoc=ExpTest_all_trials_InLoc(rand_parameter,:);
    loc_allCondsInOrder=loc_allCondsInOrder(rand_parameter,:);
    CondInOrder=reshape(repmat([1 2 3 4],pairs_per_cond,1),pairs_per_cond*4,1);
    CondInOrder=CondInOrder(rand_parameter);
    %check that a target/distractor do not repeat within two trials
    for r=2:length(ExpTest_all_trials_InLoc)
        if any(ismember(ExpTest_all_trials_InLoc(r,:),ExpTest_all_trials_InLoc(r-1,:)))
            check=1;
            break;
        end
        if r>2 && check==0 %need to check the 2-back
            if any(ismember(ExpTest_all_trials_InLoc(r,:),ExpTest_all_trials_InLoc(r-2,:)))
                check=1;
                break;
            end
        end
        
    end
    
    if mod(attempts,100000)==0
        fprintf(sprintf('attempt Exp part 1 %d\n',attempts));
    end
    if attempts > 2000000
        check=0;
        fprintf('didn''t find a good sequence for Exp test part 1, run script again\n');
    end
end
save([subj_dir '/Explicit_trial_sequencesPart1.mat'],'ExpTest_all_trials_InLoc','ExpTest_cues','ExpTest_dist','loc_allConds','loc_allCondsInOrder','CondInOrder');

%% explicit memory test - PART 2: testing for B-C v and nv
%these are the distractors, different randomization than the one for the
%rogue
smallDist=[[5:(pairs_per_cond/2) 1:4]' [7:(pairs_per_cond/2) 1:6]'];
bigDist=smallDist+(pairs_per_cond/2);
AllDist=[smallDist;bigDist];

All_stim(:,:,1)=stim_v;
All_stim(:,:,2)=stim_nv;
num_conds=size(All_stim,3);

ExpTest_cues=squeeze(All_stim(:,1,:)); %all A items
ExpTest_dist=zeros(pairs_per_cond,3,num_conds); %each row is a trial, column 1: B target. column 2: dist. coumn 3: dist. Each z level is a condition

for cond=1:num_conds
    ExpTest_dist(:,1,cond)=All_stim(:,2,cond);%all B items
end

for cond=1:num_conds
    ExpTest_dist(:,2:3,cond)=[All_stim(AllDist(:,1),2,cond) All_stim(AllDist(:,2),2,cond)];%all C items
end

%build the locations matrix:
loc=perms([1,2,3]);
all_locs=repmat(loc,3,1);%to have 18 - for big and small.
loc_allConds=zeros(pairs_per_cond,3,num_conds);
for cond=1:num_conds %each condition
    loc_allConds(:,:,cond)=[all_locs(randperm(size(all_locs,1)),:);all_locs(randperm(size(all_locs,1)),:)];
end

%build the stimuli in according to their locations
for cond=1:num_conds
    for pair=1:pairs_per_cond
    ExpTest_dist_inLoc(pair,:,cond)=ExpTest_dist(pair,loc_allConds(pair,:,cond),cond);
    end
end

%construct the order of the pairs - make sure that there's a gap between identical items;
%first three columns are the distractors, the cue is the forth column.
check=1;
attempts=0;
while check
    check=0;
    attempts=attempts+1;
    ExpTest_all_trials_InLoc=[];
    loc_allCondsInOrder=[];
    for cond=1:num_conds
        ExpTest_all_trials_InLoc=[ExpTest_all_trials_InLoc; [ExpTest_dist_inLoc(:,:,cond) ExpTest_cues(:,cond)]];
        loc_allCondsInOrder=[loc_allCondsInOrder;loc_allConds(:,:,cond)];
    end
    %shuffle:
    rand_parameter=randperm(length(ExpTest_all_trials_InLoc))';
    ExpTest_all_trials_InLoc=ExpTest_all_trials_InLoc(rand_parameter,:);
    loc_allCondsInOrder=loc_allCondsInOrder(rand_parameter,:);
    CondInOrder=reshape(repmat((1:num_conds),pairs_per_cond,1),pairs_per_cond*num_conds,1);
    CondInOrder=CondInOrder(rand_parameter);
    for r=2:length(ExpTest_all_trials_InLoc)
        if any(ismember(ExpTest_all_trials_InLoc(r,:),ExpTest_all_trials_InLoc(r-1,:)))
            check=1;
            break;
        end
        if r>2 && check==0 %need to check the 2-back
            if any(ismember(ExpTest_all_trials_InLoc(r,:),ExpTest_all_trials_InLoc(r-2,:)))
                check=1;
                break;
            end
        end
    end
    
    if mod(attempts,100000)==0
        fprintf(sprintf('attempt Exp part 2 %d\n',attempts));
    end
    if attempts > 1000000
        check=0;
        fprintf('didn''t find a good sequence for Exp test part 2, run script again\n');
    end
    
end
save([subj_dir '/Explicit_trial_sequencesPart2.mat'],'ExpTest_all_trials_InLoc','ExpTest_cues','ExpTest_dist','loc_allConds','loc_allCondsInOrder','CondInOrder');

end

 
%% subfunctions




        

