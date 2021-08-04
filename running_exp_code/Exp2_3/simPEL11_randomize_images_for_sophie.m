function simPEL11_randomize_images_for_sophie(subj_num,subj_id)

%% initial definitions
rand('state',sum(100.*clock));%solve the randperm problem in matlab from being consistent

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
categ_names={'small','big'}; %I called big/small "categories"
num_categ=length(categ_names);

foils_per_categ=[pairs_per_cond pairs_per_cond]./2;%number of foils during item test in each category: small, big - since only half of the items will be presented, number of foils should correspond to that.
num_Ssmall=pairs_per_cond/num_categ; %number of single items small in each v/nv condition - make sure that pairs_per_cond/num_categ dvide to integers
num_Sbig=pairs_per_cond/num_categ; %number of single items big in each v/nv condition - make sure that pairs_per_cond/num_categ dvide to integers
%single_per_categ=[num_Ssmall*num_conds num_Sbig*num_conds]; %total single items
items_per_cond_per_categ=pairs_per_cond*num_conds+foils_per_categ;%big/small items that will be sllocated to either foils/initial pairs in the study
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

%% prepare the randomizations of the predictive pairs of objects %%
% setting the randomization of objects' sizes
%the response wanted - -1 means smaller 1 means bigger
%there will be big/small objects in the experiment

%there are 36 orders in total: first 18 will be big-small, last 18 will be
%small-big
size_mat=[repmat([2,1],pairs_per_cond/2,1);repmat([1,2],pairs_per_cond/2,1)];
response_mat=[repmat([1,-1],pairs_per_cond/2,1);repmat([-1,1],pairs_per_cond/2,1)];

%now, for each condition, choose x items of each size, and assign them to
%pairs randomly.

%stim_alocation 9randomize big and small items):
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
%insert the wanted response into the sequence (1:bigger, -1: smaller) -
%we'll use that to determine accuracy
stim_v=stim_v.*response_mat;
stim_nv=stim_nv.*response_mat;

%foils for test - take half big and half small
stim_foils=[];
for c=1:num_categ
    stim_foils=[stim_foils; stim_rand((num_conds*pairs_per_cond+1):items_per_cond_per_categ(c),c)+9000];
end

%randomize the order:
stim_foils=stim_foils(randperm(length(stim_foils)));

%% grab items for all single items:2 conditions,in each there are small and big items

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
%                                             single violation=5000  single no-violation=6000 single violation lures=7000  single no-violation lures=8000 foil=9000.
%Then the stimul: 1-x in each size.

%save it - need the same randomization for all parts
save([subj_dir '/items_randomization'],'stim_v','stim_nv','stim_single_v','stim_single_nv','lures_single_v','lures_single_nv','stim_foils');

%% construct trials sequences for initial learning sessions

%construct pair sequence for initial learning, day 1 (v,nv):
pairs_temp=[ones(1,pairs_per_cond)*1000 ones(1,pairs_per_cond)*2000];
pairs_temp=pairs_temp+[1:pairs_per_cond 1:pairs_per_cond];

pairs_seq_day1=zeros(reps_init_learning,length(pairs_temp));

%now run a loop until the sequence is good:
counter=0; %to count how long it took, good for debugging
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
            %repetitions (cannot happen within each cycle of repetition,
            %since each pair appears once.
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
    %not sure why this is here, it doesn't make sense.. doesn't do anything...:
    stim_in_pair_seq_day1(stim_in_pair_seq_day1==0)=pairs_per_cond;%correct the modulus
    
    if counter==100000
        break
    end
    
end %ends the while loop that create the sequence for the initial learning day 1

%% generate onset times:
%that is needed for the randomization of the pairs later on:
rand_pairs_isi=randperm(pairs_per_cond);
%two rows, dividing the pairs that will appear in the two sections of the reminder by similar isi combo, to equate the length
%this is unnecessary w/o jittering, but I never removed it.
rand_pairs_rem=[rand_pairs_isi(1:2:end) rand_pairs_isi(2:2:end)];

% set trial duration and ITI: I had it set bc in previous iterations that
% included jittter for fMRI, but no longer needed
% was jittered
run_delay = 2;
run_decay = 2;
stim_dur = 1.5;
isi=0.5;
isi_day1=ones(size(stim_seq_day1))*isi;
trial_lengths=stim_dur + isi_day1;
day1_init_onsets=run_delay + cumsum(trial_lengths,2) - (stim_dur+isi_day1);


%% day2: %%%%%
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
stim_seq_reminder=zeros(runs_rem,seq_length/runs_rem*2); %sequence of stim
stim_loc_seq_reminder=zeros(size(stim_seq_reminder)); %the location within a pair (1/2)
stim_in_pair_seq_reminder=zeros(size(stim_seq_reminder)); %which pair number it is

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

%% generate onsets reminder sessions
%as above - used to  be jittered, now not, unnecessary
%parameters of the simulation
run_delay = 2;
run_decay = 2;

%% put in the jitter:

%no jitter:
isi_reminder=ones(size(stim_seq_reminder))*isi;

trial_lengths=stim_dur + isi_reminder;
reminder_onsets=run_delay + cumsum(trial_lengths,2) - (stim_dur+isi_reminder);

%% construct trials sequences for violation block

runs_violation=4; %divide to 4 blocks, each will have 9 pairs in each condition, reminder and v/nv
pairs_per_run=pairs_per_cond/runs_violation;
check_length=1;
attempts_length=0;

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
    
    %violation: will be inseated instead of second items.
    %no violation: will be insearted after the nv pair
    stim_seq_vnv=zeros(runs_violation,size(pairs_seq_vnv,2)*2+pairs_per_cond/runs_violation); %in all conditions there are 2 items per pair, and add pairs_per_cond/runs_violation more for the nv items that are added to nv pairs
    stim_loc_seq_vnv=zeros(size(stim_seq_vnv));
    stim_in_pair_seq_vnv=zeros(size(stim_seq_vnv));
    
    %%%% for the nv single items - they are not insearted instead of the B item,
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
    
    %% generate onsets 
    isi_vnv_block=ones(size(stim_seq_vnv))*isi;
    trial_lengths=stim_dur + isi_vnv_block;
    vnv_onsets=run_delay + cumsum(trial_lengths,2) - (stim_dur+isi_vnv_block);
    
    
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
end

% when running the script, check that this is okay.
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


%%
%used to be jittered, not any more, but code is as if it is jittered... :)
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

%% explicit memory test - PART 1: testing for A-C (v), B-C (nv)
%Set the distractors for the explicit memory test. 
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

%% explicit memory test - PART 2: testing for A-B v and nv
%these are the distractors:
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



