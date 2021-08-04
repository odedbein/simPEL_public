%% this script is not meant to run as a function - but by executing sections

%% I only reported accuracy baseed on this script, all else is in R.
%% But, has some ploting of all kinds of things if anyone is interested
%set stats parameters:
showtable=1;
alpha=0.05;
%set initial parameters:
items=['A','B'];
prev_chunk_size=numel(items);%in this version - pairs, in previous triads - would be 3.
num_reps_init=12;

plot_single_sub=0;
logRT=1;
num_cond=2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% all items regardless of AB memory:
%% day1_init:
subjects_numbers=[50:66 68:82];

good_subj=[1:2 4:7 9:15 17:21 23:32]; 


header.day1_init={};
for item=1:prev_chunk_size
    for rep=1:num_reps_init
    header.day1_init{1,(item-1)*num_reps_init+rep}=sprintf('violation%s%d',items(item),rep);
    header.day1_init{1,prev_chunk_size*num_reps_init+(item-1)*num_reps_init+rep}=sprintf('no violation%s%d',items(item),rep);
    end
end

header.day1_init(1,2:end+1)=header.day1_init(1,1:end);
header.day1_init{1,1}='subjects';

if logRT
    RT=simPEL11_All_subs_learning_byItemAllReps.logRT.AllItems.day1_init(good_subj,:);
else
    RT=simPEL11_All_subs_learning_byItemAllReps.RT.AllItems.day1_init(good_subj,:);
end
RTStruct.day1_init=header.day1_init;
RTStruct.day1_init(2:length(good_subj)+1,1)=num2cell(subjects_numbers(good_subj)');
RTStruct.day1_init(2:length(good_subj)+1,2:end)=num2cell(RT);

subs_num_array={};
for i=1:length(good_subj)
    subs_num_array{i}=num2str(subjects_numbers(i));
end
cmap = jet(length(good_subj));
if plot_single_sub
    for i=1:num_reps_init:size(RT,2)
        curr_data=RT(:,i:i+num_reps_init-1)';
        figure;
        for subj=1:length(good_subj)
            line(1:num_reps_init, curr_data(:,subj), 'Color', cmap(subj, :));
        end
        set(gca,'xtick',1:num_reps_init);
        Xticks=RTStruct.day1_init(1,i+1:i+num_reps_init);
        set(gca,'XTickLabel',Xticks);
        legend(subs_num_array,'Fontsize',12);
    end
end


%accuracy structure:
AccStruct.day1_init=RTStruct.day1_init;
AccStruct.day1_init(2:length(good_subj)+1,2:end)=num2cell(simPEL11_All_subs_learning_byItemAllReps.accuracy_rate.AllItems.day1_init(good_subj,:));

%calculate the average accuracy for all relevant items (B) for each
%repetition:
AccResp=simPEL11_All_subs_learning_byItemAllReps.accuracy_rate.AllItems.day1_init(good_subj,:);
AvAccResp=[];
for rep=1:num_reps_init
    AvAccResp=[AvAccResp mean(AccResp(:,[num_reps_init+rep,num_reps_init*3+rep]),2)];%num_reps_init - jump the A items in the violation condition, num_reps_init*2 - jump the A items in the nv condition
end

%reported accuracy:
meanAvAccResp_perSubj=mean(AvAccResp,2);
fprintf('accuracy during initial learning: M = %.2f, SD: %.2f\n',mean(meanAvAccResp_perSubj),std(meanAvAccResp_perSubj));


meanRT=(nanmean(RT));
meanRTABC_day1_init=reshape(meanRT,num_reps_init*length(items),numel(meanRT)/(num_reps_init*prev_chunk_size))';
meanRTABCStruct={};
meanRTABCStruct.day1_init{1,1}={'condition'};
for item=1:prev_chunk_size
    for rep=1:num_reps_init
    meanRTABCStruct.day1_init{1,1+(item-1)*num_reps_init+rep}=sprintf('%s%d',items(item),rep);
    end
end

meanRTABCStruct.day1_init(2,1)={'violation'};
meanRTABCStruct.day1_init(3,1)={'no violation'};
meanRTABCStruct.day1_init(2:3,2:end)=num2cell(meanRTABC_day1_init);

%% some ploting - that will be done only with good items. B items in each condition, compare between conditions - so SEM should be between conditions:

if logRT
    RT=simPEL11_All_subs_learning_byItemAllReps.logRT.AllItems.day1_init(good_subj,:);
else
    RT=simPEL11_All_subs_learning_byItemAllReps.RT.AllItems.day1_init(good_subj,:);
end

%plot B items in each condition
curr_data=[];
for c=1:num_cond
    curr_data=[curr_data RT(:,(((c-1)*num_reps_init*prev_chunk_size+num_reps_init+1):((c-1)*num_reps_init*prev_chunk_size+num_reps_init*2)))];
end
meanRT=nanmean(curr_data)';
meanRT=reshape(meanRT,num_reps_init,num_cond)';
SEM=[];
for i=1:num_reps_init
    subjAv=nanmean(curr_data(:,[i,num_reps_init+i]),2);
    withinEr=curr_data(:,[i,num_reps_init+i])-repmat(subjAv,1,num_cond);
    SEM=[SEM (nanstd(withinEr)/sqrt(size(curr_data,1)))'];
end

figure;
hold on
for i=1:num_cond %violation, no violation
    errorbar(1:num_reps_init,meanRT(i,:),SEM(i,:));
end
set(gca,'xtick',1:num_reps_init);
legend('violaiton','no violation');
title('B items All items','Fontsize',16)
if logRT
    ylabel('logRT','Fontsize',16);
else
    ylabel('RT(ms)','Fontsize',16);
end
xlabel('repetition','Fontsize',16);
hold off

%% stats:
%for the anova (and currently subsequent stats)

if logRT
    RT=simPEL11_All_subs_learning_byItemAllReps.logRT.AllItems.day1_init(good_subj,:);
else
    RT=simPEL11_All_subs_learning_byItemAllReps.RT.AllItems.day1_init(good_subj,:);
end

%get B items in each condition
curr_data=[];
for c=1:num_cond
    curr_data=[curr_data RT(:,(((c-1)*num_reps_init*prev_chunk_size+num_reps_init+1):((c-1)*num_reps_init*prev_chunk_size+num_reps_init*2)))];
end

n=size(curr_data,1);
Y=reshape(curr_data,size(curr_data,1)*size(curr_data,2),1);
S=repmat([1:n]',num_reps_init*2,1);
F1=[ones(n*num_reps_init,1);ones(n*num_reps_init,1)*2];%v/nv
F2=reshape(repmat(1:num_reps_init,n,2),size(Y));%repetitions
fprintf('ANOVA initial learning B items RT: F1:v/nv, F2: repetitions 1-7\n');
X=[Y,F1,F2,S];
RMAOV2_mod(X,alpha,showtable);

%t-tests for each repetition:
for irep=1:num_reps_init
    col1=irep;
    col2=irep+num_reps_init;
    fprintf('B items repetition %d: v vs. nv ttest: \n',irep);
    [h,p,ci,stats] = ttest(curr_data(:,col1),curr_data(:,col2));
    p
end

%t-tests first vs. last repetition:
col1=1;
col2=num_reps_init;
fprintf('violation: rep 1-12 \n');
[h,p,ci,stats] = ttest(curr_data(:,col1),curr_data(:,col2));
p

col1=num_reps_init+1;
col2=num_reps_init*2;
fprintf('no-violation: rep 1-12  \n');
[h,p,ci,stats] = ttest(curr_data(:,col1),curr_data(:,col2));
p

%% plot A items in each condition
if logRT
    RT=simPEL11_All_subs_learning_byItemAllReps.logRT.AllItems.day1_init(good_subj,:);
else
    RT=simPEL11_All_subs_learning_byItemAllReps.RT.AllItems.day1_init(good_subj,:);
end

curr_data=[];
for c=1:num_cond
    curr_data=[curr_data RT(:,(((c-1)*num_reps_init*length(items)+1):((c-1)*num_reps_init*length(items)+num_reps_init)))];
end

meanRT=nanmean(curr_data)';
meanRT=reshape(meanRT,num_reps_init,num_cond)';
SEM=[];
for i=1:num_reps_init
    subjAv=nanmean(curr_data(:,[i,num_reps_init+i]),2);
    withinEr=curr_data(:,[i,num_reps_init+i])-repmat(subjAv,1,num_cond);
    SEM=[SEM (nanstd(withinEr)/sqrt(size(curr_data,1)))'];
end

figure;
hold on
for i=1:num_cond %violation, no violation
    errorbar(1:num_reps_init,meanRT(i,:),SEM(i,:));
end
set(gca,'xtick',1:num_reps_init);
legend('violaiton','no violation');
title('A items','Fontsize',16)
if logRT
    ylabel('logRT','Fontsize',16);
else
    ylabel('RT(ms)','Fontsize',16);
end
xlabel('repetition','Fontsize',16);
hold off

%% stats

if logRT
    RT=simPEL11_All_subs_learning_byItemAllReps.logRT.AllItems.day1_init(good_subj,:);
else
    RT=simPEL11_All_subs_learning_byItemAllReps.RT.AllItems.day1_init(good_subj,:);
end

%get A items in each condition
curr_data=[];
for c=1:num_cond
    curr_data=[curr_data RT(:,(((c-1)*num_reps_init*length(items)+1):((c-1)*num_reps_init*length(items)+num_reps_init)))];
end


n=size(curr_data,1);
Y=reshape(curr_data,size(curr_data,1)*size(curr_data,2),1);
S=repmat([1:n]',num_reps_init*2,1);
F1=[ones(n*num_reps_init,1);ones(n*num_reps_init,1)*2];%v/nv
F2=reshape(repmat(1:num_reps_init,n,2),size(Y));%repetitions
fprintf('ANOVA initial learning A items RT: F1:v/nv, F2: repetitions 1-7\n');
X=[Y,F1,F2,S];
RMAOV2_mod(X,alpha,showtable);

%t-tests for each repetition:
for irep=1:num_reps_init
    col1=irep;
    col2=irep+num_reps_init;
    fprintf('A items repetition %d: v vs. nv ttest: \n',irep);
    [h,p,ci,stats] = ttest(curr_data(:,col1),curr_data(:,col2));
    p
end

%t-tests first vs. last repetition:
col1=1;
col2=num_reps_init;
fprintf('violation: rep 1-7 \n');
[h,p,ci,stats] = ttest(curr_data(:,col1),curr_data(:,col2));
p

col1=num_reps_init+1;
col2=num_reps_init*2;
fprintf('no-violation: rep 1-7  \n');
[h,p,ci,stats] = ttest(curr_data(:,col1),curr_data(:,col2));
p

%% compare A and B items, collapse on conditions - to see learning of B items. should be similar in both conditions, so first just plot collapsed, then see if no difference.

if logRT
    curr_data=simPEL11_All_subs_learning_byItemAllRepsCollapsVnv.logRT.AllItems.day1_init(good_subj,:);
else
    curr_data=simPEL11_All_subs_learning_byItemAllRepsCollapsVnv.RT.AllItems.day1_init(good_subj,:);
end

meanRT=nanmean(curr_data)';
meanRT=reshape(meanRT,num_reps_init,num_cond)'; %two rows, A and B
SEM=[];
for i=1:num_reps_init
    subjAv=nanmean(curr_data(:,[i,num_reps_init+i]),2);
    withinEr=curr_data(:,[i,num_reps_init+i])-repmat(subjAv,1,num_cond);
    SEM=[SEM (nanstd(withinEr)/sqrt(size(curr_data,1)))'];
end

colors=[...
        0.5 0.5 0.5
        0   0   0
        ];
    
figure;
hold on
for i=1:prev_chunk_size %A,B
    errorbar(1:num_reps_init,meanRT(i,:),SEM(i,:),'Color',colors(i,:));
end
set(gca,'xtick',1:num_reps_init);
legend('A items','B items');
title('learning day1 collapse conditions','Fontsize',16)
if logRT
    ylabel('logRT','Fontsize',16);
else
    ylabel('RT(ms)','Fontsize',16);
end
xlabel('repetition','Fontsize',16);
hold off

%% stats
 
if logRT
    curr_data=simPEL11_All_subs_learning_byItemAllRepsCollapsVnv.logRT.AllItems.day1_init(good_subj,:);
else
    curr_data=simPEL11_All_subs_learning_byItemAllRepsCollapsVnv.RT.AllItems.day1_init(good_subj,:);
end

n=size(curr_data,1);
Y=reshape(curr_data,size(curr_data,1)*size(curr_data,2),1);
S=repmat([1:n]',num_reps_init*2,1);
F1=[ones(n*num_reps_init,1);ones(n*num_reps_init,1)*2];%A/B items
F2=reshape(repmat(1:num_reps_init,n,2),size(Y));%repetitions
fprintf('ANOVA initial learning day1 RT: F1:A/B item, F2: repetitions 1-12\n');
X=[Y,F1,F2,S];
RMAOV2_mod(X,alpha,showtable);

%t-tests for each repetition:
for irep=1:num_reps_init
    col1=irep;
    col2=irep+num_reps_init;
    fprintf('A vs. B items repetition %d: v vs. nv ttest: \n',irep);
    [h,p,ci,stats] = ttest(curr_data(:,col1),curr_data(:,col2));
    p
end

%% reminder - data structures
header.rem={'subjects','violationA','violationB','no-violationA','no-violationB'};

%good_subj=1:length(subjects_numbers);
if logRT
    RT=simPEL11_All_subs_learning_byItemAllReps.logRT.AllItems.rem(:,:);
else
    RT=simPEL11_All_subs_learning_byItemAllReps.RT.AllItems.rem(:,:);
end
RTStruct.rem=header.rem;
RTStruct.rem(2:length(subjects_numbers)+1,1)=num2cell(subjects_numbers');
RTStruct.rem(2:length(subjects_numbers)+1,2:end)=num2cell(RT);
%accuracy structure:
AccStruct.rem=RTStruct.rem;
AccStruct.rem(2:length(subjects_numbers)+1,2:end)=num2cell(simPEL11_All_subs_learning_byItemAllReps.accuracy_rate.AllItems.rem(:,:));

%calculate the average accuracy for all relevant items (B)
AvAccResp=simPEL11_All_subs_learning_byItemAllReps.accuracy_rate.AllItems.rem(good_subj,[2,4]);
%reported accuracy:
meanAvAccResp_perSubj=mean(AvAccResp,2);
fprintf('accuracy during reminder: M = %.2f, SD: %.2f\n',mean(meanAvAccResp_perSubj),std(meanAvAccResp_perSubj));


meanRT=nanmean(RT);
meanRTABC_rem=[meanRT(1:length(items));meanRT(length(items)+1:length(items)*2);];
meanRTABCStruct.rem={'condition','A','B'};
meanRTABCStruct.rem(2,1)={'violation'};
meanRTABCStruct.rem(3,1)={'no-violation'};
meanRTABCStruct.rem(2:3,2:end)=num2cell(meanRTABC_rem);

%% ploting - A and B in each condition
plotSubj=0;
if logRT
    RT=simPEL11_All_subs_learning_byItemAllReps.logRT.AllItems.rem(good_subj,:);
else
    RT=simPEL11_All_subs_learning_byItemAllReps.RT.AllItems.rem(good_subj,:);
end

%plot A and B items in each condition - group by condition:
figure;
subs_col=get(gca,'colororder');

curr_data=RT;
num_subj=size(curr_data,1);
meanRT=nanmean(curr_data);
%compare within condition, A vs. B:
SEM=[];
for i=[1,3]
    subjAv=nanmean(curr_data(:,[i,i+1]),2);
    withinEr=curr_data(:,[i,i+1])-repmat(subjAv,1,num_cond);
    SEM=[SEM (nanstd(withinEr)/sqrt(size(curr_data,1)))'];
end


color=[
    0.5       0.5       0.5
    0         0.4470    0.7410
    0.5       0.5       0.5
    0.8500    0.3250    0.0980];
    

subplot(1,2,1);
bar(1:num_cond*prev_chunk_size,meanRT,'FaceColor',[1,1,1]);
hold on
for i=1:num_cond*prev_chunk_size%
    bar(i,meanRT(i),'FaceColor',color(i,:));
    errorbar(i,meanRT(i),SEM(i),'k');
end

if logRT
    ylabel('logRT','Fontsize',16);
    yLimits=[6,6.6];
else
    ylabel('RT(ms)','Fontsize',16);
    yLimits=[620,800];
end
title('simPEL11 reminder','Fontsize',20);
Xticks={'violation: A','violation: B','no-violation: A','no-violation: B'};
set(gca,'XTickLabel',Xticks,'FontSize', 14,'XTickLabelRotation',25);

% Plot subject data
if plotSubj
    subjDataX=[];
    subjDataY=[];
    for i=1:num_cond*prev_chunk_size
        subjDataX = [subjDataX ones(1,num_subj)*i];
        subjDataY = [subjDataY curr_data(:,i)'];
    end
    scatter(subjDataX,subjDataY, 15, 'ko', 'filled')
    % Plot lines connecting
    for ii = 1:num_subj
        plot(1:prev_chunk_size,curr_data(ii,1:prev_chunk_size), 'Color', subs_col((mod(ii,size(subs_col,1))+1),:))
        plot(prev_chunk_size+(1:prev_chunk_size),curr_data(ii,prev_chunk_size+(1:prev_chunk_size)), 'Color', subs_col((mod(ii,size(subs_col,1))+1),:))
    end
    if logRT
       yLimits=[6,7];
    else
       yLimits=[500,1050];
    end
end
ylim(yLimits);

hold off

col1=1;
col2=2;
fprintf('reminder A vs. B items violation: \n');
[h,p,ci,stats] = ttest(curr_data(:,col1),curr_data(:,col2));
p

col1=3;
col2=4;
fprintf('reminder A vs. B items n-violation: \n');
[h,p,ci,stats] = ttest(curr_data(:,col1),curr_data(:,col2));
p

%plot A and B items in each condition - group by item type:
cond_shuffle=[1,3,2,4];
curr_data=RT(:,cond_shuffle);
num_subj=size(curr_data,1);
meanRT=nanmean(curr_data);
%compare within condition, A vs. B:
SEM=[];
for i=[1,3]
    subjAv=nanmean(curr_data(:,[i,i+1]),2);
    withinEr=curr_data(:,[i,i+1])-repmat(subjAv,1,num_cond);
    SEM=[SEM (nanstd(withinEr)/sqrt(size(curr_data,1)))'];
end

color=color(cond_shuffle,:);
    
subplot(1,2,2);
bar(1:num_cond*prev_chunk_size,meanRT,'FaceColor',[1,1,1]);
hold on
for i=1:num_cond*prev_chunk_size%
    bar(i,meanRT(i),'FaceColor',color(i,:));
    errorbar(i,meanRT(i),SEM(i),'k');
end

if logRT
    ylabel('logRT','Fontsize',16);
    yLimits=[6,6.6];
else
    ylabel('RT(ms)','Fontsize',16);
    yLimits=[620,800];
end
title('simPEL11 reminder','Fontsize',20);
set(gca,'XTickLabel',Xticks(cond_shuffle),'FontSize', 14,'XTickLabelRotation',25);

% Plot subject data
if plotSubj
    subjDataX=[];
    subjDataY=[];
    for i=1:num_cond*prev_chunk_size
        subjDataX = [subjDataX ones(1,num_subj)*i];
        subjDataY = [subjDataY curr_data(:,i)'];
    end
    scatter(subjDataX,subjDataY, 15, 'ko', 'filled')
    % Plot lines connecting
    for ii = 1:num_subj
        plot(1:prev_chunk_size,curr_data(ii,1:prev_chunk_size), 'Color', subs_col((mod(ii,size(subs_col,1))+1),:))
        plot(prev_chunk_size+(1:prev_chunk_size),curr_data(ii,prev_chunk_size+(1:prev_chunk_size)), 'Color', subs_col((mod(ii,size(subs_col,1))+1),:))
    end
    if logRT
       yLimits=[6,7];
    else
       yLimits=[500,1050];
    end
end
ylim(yLimits);

hold off

col1=1;
col2=2;
fprintf('reminder A items: v vs. nv: \n');
[h,p,ci,stats] = ttest(curr_data(:,col1),curr_data(:,col2));
p

col1=3;
col2=4;
fprintf('reminder B items: v vs. nv: \n');
[h,p,ci,stats] = ttest(curr_data(:,col1),curr_data(:,col2));
p

           
%% error block - data structures
num_cond=2;
header.er_block={'subjects','violA-rem','violB-rem','no-violA-rem','no-violB-rem','violA-viol','violB-viol','no-violA-viol','no-violB-viol','viol-single','no-viol-single'};
logRT=0;
if logRT
    RT=simPEL11_All_subs_learning_byItemAllReps.logRT.AllItems.er_block(:,:);
else
    RT=simPEL11_All_subs_learning_byItemAllReps.RT.AllItems.er_block(:,:);
end
RTStruct.er_block=header.er_block;
RTStruct.er_block(2:length(subjects_numbers)+1,1)=num2cell(subjects_numbers');
RTStruct.er_block(2:length(subjects_numbers)+1,2:end)=num2cell(RT);
%accuracy structure:
AccStruct.er_block=RTStruct.er_block;
AccStruct.er_block(2:length(subjects_numbers)+1,2:end)=num2cell(simPEL11_All_subs_learning_byItemAllReps.accuracy_rate.AllItems.er_block(:,:));

%calculate the average accuracy for all relevant items (B)
AvAccResp=simPEL11_All_subs_learning_byItemAllReps.accuracy_rate.AllItems.er_block(good_subj,[2,4,8:10]);

%reported accuracy:
meanAvAccResp_perSubj=mean(AvAccResp,2);
fprintf('accuracy during violation: M = %.2f, SD: %.2f\n',mean(meanAvAccResp_perSubj),std(meanAvAccResp_perSubj));
%compare violation and no violation:
fprintf('accuracy v items: M = %.2f, SD: %.2f\n',mean(AvAccResp(:,4)),std(AvAccResp(:,4)));
fprintf('accuracy nv items: M = %.2f, SD: %.2f\n',mean(AvAccResp(:,5)),std(AvAccResp(:,5)));
fprintf('accuracy violation: v vs. nv items: \n');
[h,p,ci,stats] = ttest(AvAccResp(:,4),AvAccResp(:,5))


meanRT=mean(RT);
meanRTABC_er_block=[meanRT(1:length(items));meanRT(length(items)+1:length(items)*2);meanRT(length(items)*2+1:length(items)*3);meanRT(length(items)*3+1:length(items)*4)];
meanRTABCStruct.er_block={'condition','A','B'};
meanRTABCStruct.er_block(2,1)={'viol-rem'};
meanRTABCStruct.er_block(3,1)={'no-viol-rem'};
meanRTABCStruct.er_block(4,1)={'viol-viol'};
meanRTABCStruct.er_block(5,1)={'no-viol-viol'};
meanRTABCStruct.er_block(2:5,2:end)=num2cell(meanRTABC_er_block);


%% ploting - A and B in the reminders during the v/nv block condition
if logRT
    RT=simPEL11_All_subs_learning_byItemAllReps.logRT.AllItems.er_block(good_subj,:);
else
    RT=simPEL11_All_subs_learning_byItemAllReps.RT.AllItems.er_block(good_subj,:);
end

%plot A and B items in each condition - group by condition:
figure;
subs_col=get(gca,'colororder');

curr_data=RT(:,1:4);
num_subj=size(curr_data,1);
meanRT=nanmean(curr_data);
%compare within condition, A vs. B:
SEM=[];
for i=[1,3]
    subjAv=nanmean(curr_data(:,[i,i+1]),2);
    withinEr=curr_data(:,[i,i+1])-repmat(subjAv,1,num_cond);
    SEM=[SEM (nanstd(withinEr)/sqrt(size(curr_data,1)))'];
end


color=[
    0.5       0.5       0.5
    0         0.4470    0.7410
    0.5       0.5       0.5
    0.8500    0.3250    0.0980];
    

subplot(1,2,1);
bar(1:num_cond*prev_chunk_size,meanRT,'FaceColor',[1,1,1]);
hold on
for i=1:num_cond*prev_chunk_size%
    bar(i,meanRT(i),'FaceColor',color(i,:));
    errorbar(i,meanRT(i),SEM(i),'k');
end

if logRT
    ylabel('logRT','Fontsize',16);
    yLimits=[6,6.6];
else
    ylabel('RT(ms)','Fontsize',16);
    yLimits=[620,800];
end
title('simPEL11 v/nv block: reminder','Fontsize',20);
Xticks={'violation: A','violation: B','no-violation: A','no-violation: B'};
set(gca,'XTickLabel',Xticks,'FontSize', 14,'XTickLabelRotation',25);

% Plot subject data
if plotSubj
    subjDataX=[];
    subjDataY=[];
    for i=1:num_cond*prev_chunk_size
        subjDataX = [subjDataX ones(1,num_subj)*i];
        subjDataY = [subjDataY curr_data(:,i)'];
    end
    scatter(subjDataX,subjDataY, 15, 'ko', 'filled')
    % Plot lines connecting
    for ii = 1:num_subj
        plot(1:prev_chunk_size,curr_data(ii,1:prev_chunk_size), 'Color', subs_col((mod(ii,size(subs_col,1))+1),:))
        plot(prev_chunk_size+(1:prev_chunk_size),curr_data(ii,prev_chunk_size+(1:prev_chunk_size)), 'Color', subs_col((mod(ii,size(subs_col,1))+1),:))
    end
    if logRT
       yLimits=[6,7];
    else
       yLimits=[500,1050];
    end
end
ylim(yLimits);

hold off

col1=1;
col2=2;
fprintf('v/nv block: reminder A vs. B items violation: \n');
[h,p,ci,stats] = ttest(curr_data(:,col1),curr_data(:,col2));
p

col1=3;
col2=4;
fprintf('v/nv block: reminder A vs. B items n-violation: \n');
[h,p,ci,stats] = ttest(curr_data(:,col1),curr_data(:,col2));
p

%plot A and B items in each condition - group by item type:
cond_shuffle=[1,3,2,4];
curr_data=RT(:,cond_shuffle);
num_subj=size(curr_data,1);
meanRT=nanmean(curr_data);
%compare within condition, A vs. B:
SEM=[];
for i=[1,3]
    subjAv=nanmean(curr_data(:,[i,i+1]),2);
    withinEr=curr_data(:,[i,i+1])-repmat(subjAv,1,num_cond);
    SEM=[SEM (nanstd(withinEr)/sqrt(size(curr_data,1)))'];
end

color=color(cond_shuffle,:);
    
subplot(1,2,2);
bar(1:num_cond*prev_chunk_size,meanRT,'FaceColor',[1,1,1]);
hold on
for i=1:num_cond*prev_chunk_size%
    bar(i,meanRT(i),'FaceColor',color(i,:));
    errorbar(i,meanRT(i),SEM(i),'k');
end

if logRT
    ylabel('logRT','Fontsize',16);
    yLimits=[6,6.6];
else
    ylabel('RT(ms)','Fontsize',16);
    yLimits=[620,800];
end
title('simPEL11 v/nv block: reminder','Fontsize',20);
set(gca,'XTickLabel',Xticks(cond_shuffle),'FontSize', 14,'XTickLabelRotation',25);

% Plot subject data
if plotSubj
    subjDataX=[];
    subjDataY=[];
    for i=1:num_cond*prev_chunk_size
        subjDataX = [subjDataX ones(1,num_subj)*i];
        subjDataY = [subjDataY curr_data(:,i)'];
    end
    scatter(subjDataX,subjDataY, 15, 'ko', 'filled')
    % Plot lines connecting
    for ii = 1:num_subj
        plot(1:prev_chunk_size,curr_data(ii,1:prev_chunk_size), 'Color', subs_col((mod(ii,size(subs_col,1))+1),:))
        plot(prev_chunk_size+(1:prev_chunk_size),curr_data(ii,prev_chunk_size+(1:prev_chunk_size)), 'Color', subs_col((mod(ii,size(subs_col,1))+1),:))
    end
    if logRT
       yLimits=[6,7];
    else
       yLimits=[500,1050];
    end
end
ylim(yLimits);

hold off

col1=1;
col2=2;
fprintf('v/nv block: reminder A items: v vs. nv: \n');
[h,p,ci,stats] = ttest(curr_data(:,col1),curr_data(:,col2));
p

col1=3;
col2=4;
fprintf('v/nv block: reminder B items: v vs. nv: \n');
[h,p,ci,stats] = ttest(curr_data(:,col1),curr_data(:,col2));
p

%% now plot the violation part:
plotSubj=0;
logRT=1;
if logRT
    RT=simPEL11_All_subs_learning_byItemAllReps.logRT.AllItems.er_block(good_subj,:);
else
    RT=simPEL11_All_subs_learning_byItemAllReps.RT.AllItems.er_block(good_subj,:);
end

%plot A B and single items in each condition - group by condition:
figure;
subs_col=get(gca,'colororder');

curr_data=RT(:,[5 9 7 8 10]);
num_subj=size(curr_data,1);
meanRT=nanmean(curr_data);
%just all conditions, doesn't make much sense to compare here different
%item types:
subjAv=mean(curr_data,2);
withinEr=curr_data-repmat(subjAv,1,size(curr_data,2));
SEM=std(withinEr)/sqrt(size(curr_data,1));


color=[
    0.5       0.5       0.5
    0         0.2470    0.5410
    0.5       0.5       0.5
    0.8500    0.3250    0.0980
    0.6500    0.1250    0.0980];
    

subplot(1,2,1);
bar(1:size(curr_data,2),meanRT,'FaceColor',[1,1,1]);
hold on
for i=1:size(curr_data,2)%
    bar(i,meanRT(i),'FaceColor',color(i,:));
    errorbar(i,meanRT(i),SEM(i),'k');
end

if logRT
    ylabel('logRT','Fontsize',16);
    yLimits=[6,6.8];
else
    ylabel('RT(ms)','Fontsize',16);
    yLimits=[620,850];
end
title('simPEL11 v/nv block: violation','Fontsize',20);
Xticks={'violation: A','violation: single','no-violation: A','no-violation: B','no-violation: single'};
set(gca,'XTickLabel',Xticks,'FontSize', 14,'XTickLabelRotation',25);

% Plot subject data
if plotSubj
    subjDataX=[];
    subjDataY=[];
    for i=1:size(curr_data,2)
        subjDataX = [subjDataX ones(1,num_subj)*i];
        subjDataY = [subjDataY curr_data(:,i)'];
    end
    scatter(subjDataX,subjDataY, 15, 'ko', 'filled')
    % Plot lines connecting
    for ii = 1:num_subj
        plot(1:size(curr_data,2),curr_data(ii,:), 'Color', subs_col((mod(ii,size(subs_col,1))+1),:))
    end
    if logRT
       yLimits=[6,7];
    else
       yLimits=[400,1050];
    end
end
ylim(yLimits);

hold off

%plot single items in each condition:
cond_shuffle=9:10;
curr_data=RT(:,cond_shuffle);
num_subj=size(curr_data,1);
meanRT=nanmean(curr_data);
%compare sinalge items, v/nv:
subjAv=mean(curr_data,2);
withinEr=curr_data-repmat(subjAv,1,num_cond);
SEM=std(withinEr)/sqrt(size(curr_data,1));

color=color([2,5],:);
    
subplot(1,2,2);
bar(1:num_cond,meanRT,'FaceColor',[1,1,1]);
hold on
for i=1:num_cond%
    bar(i,meanRT(i),'FaceColor',color(i,:));
    errorbar(i,meanRT(i),SEM(i),'k');
end

if logRT
    ylabel('logRT','Fontsize',16);
    yLimits=[6.5,6.6];
else
    ylabel('RT(ms)','Fontsize',16);
    yLimits=[750,850];
end
title('simPEL11 v/nv block: violation','Fontsize',20);
set(gca,'XTickLabel',Xticks([2,5]),'FontSize', 14,'XTickLabelRotation',25);

% Plot subject data
if plotSubj
    subjDataX=[];
    subjDataY=[];
    for i=1:num_cond
        subjDataX = [subjDataX ones(1,num_subj)*i];
        subjDataY = [subjDataY curr_data(:,i)'];
    end
    scatter(subjDataX,subjDataY, 15, 'ko', 'filled')
    % Plot lines connecting
    for ii = 1:num_subj
        plot(1:num_cond,curr_data(ii,1:num_cond), 'Color', subs_col((mod(ii,size(subs_col,1))+1),:))
    end
    if logRT
       yLimits=[6.3,7];
    else
       yLimits=[500,1050];
    end
end
ylim(yLimits);

hold off

col1=1;
col2=2;
fprintf('v/nv block: single items: v vs. nv: \n');
[h,p,ci,stats] = ttest(curr_data(:,col1),curr_data(:,col2));
p


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BY AB MEMORY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% day1_init:
subjects_numbers=[50:66 68:82];

good_subj=[1:2 4:7 9:15 17:21 23:32]; 


header.day1_init={};
%good_subj=1:length(subjects_numbers); %just to create the structure - take all subjects
for item=1:prev_chunk_size
    for rep=1:num_reps_init
    header.day1_init{1,(item-1)*num_reps_init+rep}=sprintf('ABrem_violation%s%d',items(item),rep);
    header.day1_init{1,prev_chunk_size*num_reps_init+(item-1)*num_reps_init+rep}=sprintf('ABrem_no violation%s%d',items(item),rep);
    header.day1_init{1,prev_chunk_size*num_reps_init*2 + (item-1)*num_reps_init+rep}=sprintf('ABforg_violation%s%d',items(item),rep);
    header.day1_init{1,prev_chunk_size*num_reps_init*2 + prev_chunk_size*num_reps_init+(item-1)*num_reps_init+rep}=sprintf('ABforg_no violation%s%d',items(item),rep);
    end
end

header.day1_init(1,2:end+1)=header.day1_init(1,1:end);
header.day1_init{1,1}='subjects';

if logRT
    RT=simPEL11_All_subs_learning_byItemAllReps.logRT.Brem.day1_init(good_subj,:);
    RT=[RT simPEL11_All_subs_learning_byItemAllReps.logRT.Bforg.day1_init(good_subj,:)];
else
    RT=simPEL11_All_subs_learning_byItemAllReps.RT.Brem.day1_init(good_subj,:);
    RT=[RT simPEL11_All_subs_learning_byItemAllReps.RT.Bforg.day1_init(good_subj,:)];
end
RTStruct.day1_init=header.day1_init;
RTStruct.day1_init(2:length(good_subj)+1,1)=num2cell(subjects_numbers(good_subj)');
RTStruct.day1_init(2:length(good_subj)+1,2:end)=num2cell(RT);

subs_num_array={};
for i=1:length(good_subj)
    subs_num_array{i}=num2str(subjects_numbers(i));
end
cmap = jet(length(good_subj));
if plot_single_sub
    for i=1:num_reps_init:size(RT,2)
        curr_data=RT(:,i:i+num_reps_init-1)';
        figure;
        for subj=1:length(good_subj)
            line(1:num_reps_init, curr_data(:,subj), 'Color', cmap(subj, :));
        end
        set(gca,'xtick',1:num_reps_init);
        Xticks=RTStruct.day1_init(1,i+1:i+num_reps_init);
        set(gca,'XTickLabel',Xticks);
        legend(subs_num_array,'Fontsize',12);
    end
end

%% some ploting - that will be done only with good items. B items in each condition, compare between conditions - so SEM should be between conditions:
logRT=0;
if logRT
    RT=simPEL11_All_subs_learning_byItemAllReps.logRT.Brem.day1_init(good_subj,:);
    RT=[RT simPEL11_All_subs_learning_byItemAllReps.logRT.Bforg.day1_init(good_subj,:)];
else
    RT=simPEL11_All_subs_learning_byItemAllReps.RT.Brem.day1_init(good_subj,:);
    RT=[RT simPEL11_All_subs_learning_byItemAllReps.RT.Bforg.day1_init(good_subj,:)];
end

colors=[
    0    0.2470    0.5410
    0.6500    0.1250    0.0980
    0    0.4470    0.7410
    0.8500    0.3250    0.0980
    ];

%plot B items in each condition
curr_data=[];
for c=1:num_cond*2
    %this line takes just the B items:
    curr_data=[curr_data RT(:,(((c-1)*num_reps_init*prev_chunk_size+num_reps_init+1):((c-1)*num_reps_init*prev_chunk_size+num_reps_init*2)))];
end
meanRT=nanmean(curr_data)';
meanRT=reshape(meanRT,num_reps_init,num_cond*2)';
SEM=[];
for i=1:num_reps_init
    subjAv=nanmean(curr_data(:,[i,num_reps_init+i,num_reps_init*2+i,num_reps_init*3+i]),2);
    withinEr=curr_data(:,[i,num_reps_init+i,num_reps_init*2+i,num_reps_init*3+i])-repmat(subjAv,1,num_cond*2);
    SEM=[SEM (nanstd(withinEr)/sqrt(size(curr_data,1)))'];
end

figure;
hold on
for i=1:num_cond*2 %violation_ABrem, no violation_ABrem,violation_ABforg, no violation_ABforg
    errorbar(1:num_reps_init,meanRT(i,:),SEM(i,:),'Color',colors(i,:));
end
set(gca,'xtick',1:num_reps_init);
legend('Brem: violaiton','Brem: no violation','Bforg: violaiton','Bforg: no violation');
title('B items AB rem/forg','Fontsize',16)
if logRT
    ylabel('logRT','Fontsize',16);
else
    ylabel('RT(ms)','Fontsize',16);
end
xlabel('repetition','Fontsize',16);
hold off

%% plot A items in each condition

if logRT
    RT=simPEL11_All_subs_learning_byItemAllReps.logRT.Brem.day1_init(good_subj,:);
    RT=[RT simPEL11_All_subs_learning_byItemAllReps.logRT.Bforg.day1_init(good_subj,:)];
else
    RT=simPEL11_All_subs_learning_byItemAllReps.RT.Brem.day1_init(good_subj,:);
    RT=[RT simPEL11_All_subs_learning_byItemAllReps.RT.Bforg.day1_init(good_subj,:)];
end

curr_data=[];
for c=1:num_cond*2
    curr_data=[curr_data RT(:,(((c-1)*num_reps_init*length(items)+1):((c-1)*num_reps_init*length(items)+num_reps_init)))];
end

meanRT=nanmean(curr_data)';
meanRT=reshape(meanRT,num_reps_init,num_cond*2)';
SEM=[];
for i=1:num_reps_init
    subjAv=nanmean(curr_data(:,[i,num_reps_init+i,num_reps_init*2+i,num_reps_init*3+i]),2);
    withinEr=curr_data(:,[i,num_reps_init+i,num_reps_init*2+i,num_reps_init*3+i])-repmat(subjAv,1,num_cond*2);
    SEM=[SEM (nanstd(withinEr)/sqrt(size(curr_data,1)))'];
end

figure;
hold on
for i=1:num_cond*2 %violation_ABrem, no violation_ABrem,violation_ABforg, no violation_ABforg
    errorbar(1:num_reps_init,meanRT(i,:),SEM(i,:));
end
set(gca,'xtick',1:num_reps_init);
legend('Brem: violaiton','Brem: no violation','Bforg: violaiton','Bforg: no violation');
title('A items AB rem/forg','Fontsize',16)
if logRT
    ylabel('logRT','Fontsize',16);
else
    ylabel('RT(ms)','Fontsize',16);
end
xlabel('repetition','Fontsize',16);
hold off


%% compare A and B items, collapse on conditions - to see learning of B items. should be similar in both conditions, so first just plot collapsed, then see if no difference.

if logRT
    curr_data=simPEL11_All_subs_learning_byItemAllRepsCollapsVnv.logRT.Brem.day1_init(good_subj,:);
    curr_data=[curr_data simPEL11_All_subs_learning_byItemAllRepsCollapsVnv.logRT.Bforg.day1_init(good_subj,:)];
else
    curr_data=simPEL11_All_subs_learning_byItemAllRepsCollapsVnv.RT.Brem.day1_init(good_subj,:);
    curr_data=[curr_data simPEL11_All_subs_learning_byItemAllRepsCollapsVnv.RT.Bforg.day1_init(good_subj,:)];
end

meanRT=nanmean(curr_data)';
meanRT=reshape(meanRT,num_reps_init,num_cond*2)'; %4 rows, Brem: A and B Bforg: A and B
SEM=[];
for i=1:num_reps_init
    subjAv=nanmean(curr_data(:,[i,num_reps_init+i,num_reps_init*2+i,num_reps_init*3+i]),2);
    withinEr=curr_data(:,[i,num_reps_init+i,num_reps_init*2+i,num_reps_init*3+i])-repmat(subjAv,1,num_cond*2);
    SEM=[SEM (nanstd(withinEr)/sqrt(size(curr_data,1)))'];
end


colors=[
    0    0.2470    0.5410
    0.6500    0.1250    0.0980
    0    0.4470    0.7410
    0.8500    0.3250    0.0980
    ];

figure;
hold on
for i=1:prev_chunk_size*2 %A,B
    errorbar(1:num_reps_init,meanRT(i,:),SEM(i,:),'Color',colors(i,:));
end
set(gca,'xtick',1:num_reps_init);
legend('ABrem: A items','ABrem: B items','ABforg: A items','ABforg: B items');
title('learning day1 collapse conditions, ABrem/forg','Fontsize',16)
if logRT
    ylabel('logRT','Fontsize',16);
else
    ylabel('RT(ms)','Fontsize',16);
end
xlabel('repetition','Fontsize',16);
hold off

%% REMINDER ploting - A and B in each condition
plotSubj=0;
if logRT
    RT=simPEL11_All_subs_learning_byItemAllReps.logRT.Brem.rem(good_subj,:);
    RT=[RT simPEL11_All_subs_learning_byItemAllReps.logRT.Bforg.rem(good_subj,:)];
else
    RT=simPEL11_All_subs_learning_byItemAllReps.RT.Brem.rem(good_subj,:);
    RT=[RT simPEL11_All_subs_learning_byItemAllReps.RT.Bforg.rem(good_subj,:)];
end

Xticks={'Brem: violation: A','Brem: violation: B','Brem: no-violation: A','Brem: no-violation: B',...
        'Bforg: violation: A','Bforg: violation: B','Bforg: no-violation: A','Bforg: no-violation: B'};

color=[
    0.2       0.2       0.2
    0    0.2470    0.5410
    0.2       0.2       0.2
    0.6500    0.1250    0.0980
    0.5       0.5       0.5
    0         0.4470    0.7410
    0.5       0.5       0.5
    0.8500    0.3250    0.0980];


%plot A and B items in each condition - group by condition:
shuffle_cond=[1 5 3 7 2 6 4 8];
Xticks = Xticks(shuffle_cond);
color=color(shuffle_cond,:);

figure;
subs_col=get(gca,'colororder');

curr_data=RT(:,shuffle_cond);
num_subj=size(curr_data,1);
meanRT=nanmean(curr_data);
%compare within condition, A vs. B:
SEM=[];
for i=[1,5]
    subjAv=nanmean(curr_data(:,i:(i+3)),2);
    withinEr=curr_data(:,i:(i+3))-repmat(subjAv,1,num_cond*2);
    SEM=[SEM (nanstd(withinEr)/sqrt(size(curr_data,1)))'];
end

bar(1:length(meanRT),zeros(size(meanRT)),'FaceColor','none');
hold on
for i=1:length(meanRT)
    bar(i,meanRT(i),'FaceColor',color(i,:));
    errorbar(i,meanRT(i),SEM(i),'k');
end

if logRT
    ylabel('logRT','Fontsize',16);
    yLimits=[6.3,6.6];
else
    ylabel('RT(ms)','Fontsize',16);
    yLimits=[580,760];
end
title('simPEL11 reminder, by ABrem/forg','Fontsize',20);

set(gca,'XTickLabel',Xticks,'FontSize', 14,'XTickLabelRotation',25);

% Plot subject data
if plotSubj
    subjDataX=[];
    subjDataY=[];
    for i=1:num_cond*prev_chunk_size
        subjDataX = [subjDataX ones(1,num_subj)*i];
        subjDataY = [subjDataY curr_data(:,i)'];
    end
    scatter(subjDataX,subjDataY, 15, 'ko', 'filled')
    % Plot lines connecting
    for ii = 1:num_subj
        plot(1:prev_chunk_size,curr_data(ii,1:prev_chunk_size), 'Color', subs_col((mod(ii,size(subs_col,1))+1),:))
        plot(prev_chunk_size+(1:prev_chunk_size),curr_data(ii,prev_chunk_size+(1:prev_chunk_size)), 'Color', subs_col((mod(ii,size(subs_col,1))+1),:))
    end
    if logRT
       yLimits=[6,7];
    else
       yLimits=[500,1050];
    end
end
ylim(yLimits);

hold off

col1=1;
col2=2;
fprintf('reminder, violation A: rem-forg: \n');
[h,p,ci,stats] = ttest(curr_data(:,col1),curr_data(:,col2));
p

col1=3;
col2=4;
fprintf('reminder, no-violation A: rem-forg: \n');
[h,p,ci,stats] = ttest(curr_data(:,col1),curr_data(:,col2));
p

%remove one ss with nans:
stats_data=curr_data([1:17 19:end],1:4);
num_cond=size(stats_data,2)/2;
n=size(stats_data,1);
Y=reshape(stats_data,[numel(stats_data),1]);
S=repmat((1:n)',num_cond*2,1);
F1=[ones(n*num_cond,1);ones(n*num_cond,1)*2];%violation/no-violation 
F2=reshape(repmat(1:num_cond,n,2),size(Y));%Brem/Bforg
fprintf('ANOVA single items: F1: violation/no-violation, F2: AB memory: rem/forg \n');
X=[Y,F1,F2,S];
RMAOV2_mod(X,alpha,showtable);


col1=5;
col2=6;
fprintf('reminder, violation B: rem-forg: \n');
[h,p,ci,stats] = ttest(curr_data(:,col1),curr_data(:,col2));
p

col1=7;
col2=8;
fprintf('reminder, no-violation B: rem-forg: \n');
[h,p,ci,stats] = ttest(curr_data(:,col1),curr_data(:,col2));
p
   

%% ploting - A and B in the reminders DURING the v/nv block condition
if logRT
    RT=simPEL11_All_subs_learning_byItemAllReps.logRT.Brem.er_block(good_subj,:);
    RT=[RT simPEL11_All_subs_learning_byItemAllReps.logRT.Bforg.er_block(good_subj,:)];
else
    RT=simPEL11_All_subs_learning_byItemAllReps.RT.Brem.er_block(good_subj,:);
    RT=[RT simPEL11_All_subs_learning_byItemAllReps.RT.Bforg.er_block(good_subj,:)];
end

Xticks={'Brem: violation: A','Brem: violation: B','Brem: no-violation: A','Brem: no-violation: B',...
        'Bforg: violation: A','Bforg: violation: B','Bforg: no-violation: A','Bforg: no-violation: B'};

color=[
    0.2       0.2       0.2
    0    0.2470    0.5410
    0.2       0.2       0.2
    0.6500    0.1250    0.0980
    0.5       0.5       0.5
    0         0.4470    0.7410
    0.5       0.5       0.5
    0.8500    0.3250    0.0980];

%take only reminder:
select_cond=[1:4 11:14];
%plot A and B items in each condition - group by condition:
shuffle_cond=[1 5 3 7 2 6 4 8];
Xticks = Xticks(shuffle_cond);
color=color(shuffle_cond,:);
%plot A and B items in each condition - group by condition:
figure;
subs_col=get(gca,'colororder');

curr_data=RT(:,select_cond(shuffle_cond));
%to remove the subject without any forgotten alreafy here:
curr_data=RT([1:17 19:end],select_cond(shuffle_cond));
num_subj=size(curr_data,1);
meanRT=nanmean(curr_data);
%compare within condition, A vs. B:
for i=[1,5]
    subjAv=nanmean(curr_data(:,i:(i+3)),2);
    withinEr=curr_data(:,i:(i+3))-repmat(subjAv,1,num_cond*2);
    SEM=[SEM (nanstd(withinEr)/sqrt(size(curr_data,1)))'];
end

bar(1:length(meanRT),zeros(size(meanRT)),'FaceColor','none');
hold on
for i=1:length(meanRT)
    bar(i,meanRT(i),'FaceColor',color(i,:));
    errorbar(i,meanRT(i),SEM(i),'k');
end

if logRT
    ylabel('logRT','Fontsize',16);
    yLimits=[6.3,6.6];
else
    ylabel('RT(ms)','Fontsize',16);
    yLimits=[580,760];
end
title('simPEL11 reminder, by ABrem/forg','Fontsize',20);

set(gca,'XTickLabel',Xticks,'FontSize', 14,'XTickLabelRotation',25);

% Plot subject data
if plotSubj
    subjDataX=[];
    subjDataY=[];
    for i=1:num_cond*prev_chunk_size
        subjDataX = [subjDataX ones(1,num_subj)*i];
        subjDataY = [subjDataY curr_data(:,i)'];
    end
    scatter(subjDataX,subjDataY, 15, 'ko', 'filled')
    % Plot lines connecting
    for ii = 1:num_subj
        plot(1:prev_chunk_size,curr_data(ii,1:prev_chunk_size), 'Color', subs_col((mod(ii,size(subs_col,1))+1),:))
        plot(prev_chunk_size+(1:prev_chunk_size),curr_data(ii,prev_chunk_size+(1:prev_chunk_size)), 'Color', subs_col((mod(ii,size(subs_col,1))+1),:))
    end
    if logRT
       yLimits=[6,7];
    else
       yLimits=[500,1050];
    end
end
ylim(yLimits);

hold off

col1=1;
col2=2;
fprintf('v/nv block: reminder, violation A: rem-forg: \n');
[h,p,ci,stats] = ttest(curr_data(:,col1),curr_data(:,col2));
p

col1=3;
col2=4;
fprintf('v/nv block: reminder, no-violation A: rem-forg: \n');
[h,p,ci,stats] = ttest(curr_data(:,col1),curr_data(:,col2));
p

col1=5;
col2=6;
fprintf('v/nv block: reminder, violation B: rem-forg: \n');
[h,p,ci,stats] = ttest(curr_data(:,col1),curr_data(:,col2));
p

col1=7;
col2=8;
fprintf('v/nv block: reminder, no-violation B: rem-forg: \n');
[h,p,ci,stats] = ttest(curr_data(:,col1),curr_data(:,col2));
p

%remove one ss with nans:
stats_data=curr_data([1:17 19:end],5:8);
stats_data=curr_data(:,5:8);
num_cond=size(stats_data,2)/2;
n=size(stats_data,1);
Y=reshape(stats_data,[numel(stats_data),1]);
S=repmat((1:n)',num_cond*2,1);
F1=[ones(n*num_cond,1);ones(n*num_cond,1)*2];%violation/no-violation 
F2=reshape(repmat(1:num_cond,n,2),size(Y));%Brem/Bforg
fprintf('ANOVA Bitems: F1: violation/no-violation, F2: AB memory: rem/forg \n');
X=[Y,F1,F2,S];
RMAOV2_mod(X,alpha,showtable);

%% now plot the violation part:
plotSubj=0;
if logRT
    RT=simPEL11_All_subs_learning_byItemAllReps.logRT.Brem.er_block(good_subj,:);
    RT=[RT simPEL11_All_subs_learning_byItemAllReps.logRT.Bforg.er_block(good_subj,:)];
else
    RT=simPEL11_All_subs_learning_byItemAllReps.RT.Brem.er_block(good_subj,:);
    RT=[RT simPEL11_All_subs_learning_byItemAllReps.RT.Bforg.er_block(good_subj,:)];
end

Xticks={'Brem: violation: A','Brem: no-violation: A','Brem: no-violation: B','Brem: violation: single','Brem: no-violation: single',...
        'Bforg: violation: A','Bforg: no-violation: A','Bforg: no-violation: B','Bforg: violation: single','Bforg: no-violation: single'};

color=[
    0.2       0.2       0.2
    0.2       0.2       0.2
    0         0         0
    0         0.2470    0.5410
    0.6500    0.1250    0.0980
    
    0.5       0.5       0.5
    0.5       0.5       0.5
    0.8       0.8       0.8
    0         0.4470    0.7410
    0.8500    0.3250    0.0980];

%take only violation:
select_cond=[5 7:10 15 17:20];
%plot A and B items in each condition - group by condition:
shuffle_cond=[1 6 2 7 3 8 4 9 5 10];
Xticks = Xticks(shuffle_cond);
color=color(shuffle_cond,:);


%plot A B and single items in each condition - group by condition:
figure;
subs_col=get(gca,'colororder');

curr_data=RT(:,select_cond(shuffle_cond));
%to select participants without RT
%curr_data=RT([1:17 19:end],select_cond(shuffle_cond));
num_subj=size(curr_data,1);
meanRT=nanmean(curr_data);
%just all conditions, doesn't make much sense to compare here different
%item types:
subjAv=nanmean(curr_data,2);
withinEr=curr_data-repmat(subjAv,1,size(curr_data,2));
SEM=nanstd(withinEr)/sqrt(size(curr_data,1));
subplot(1,2,1);
bar(1:length(meanRT),zeros(size(meanRT)),'FaceColor','none');
hold on
for i=1:size(curr_data,2)%
    bar(i,meanRT(i),'FaceColor',color(i,:));
    errorbar(i,meanRT(i),SEM(i),'k');
end

if logRT
    ylabel('logRT','Fontsize',16);
    yLimits=[6.3,6.6];
else
    ylabel('RT(ms)','Fontsize',16);
    yLimits=[580,760];
end
title('simPEL11 v/nv block ABrem/forg: violation','Fontsize',20);
set(gca,'XTickLabel',Xticks,'FontSize', 14,'XTickLabelRotation',25);

% Plot subject data
if plotSubj
    subjDataX=[];
    subjDataY=[];
    for i=1:size(curr_data,2)
        subjDataX = [subjDataX ones(1,num_subj)*i];
        subjDataY = [subjDataY curr_data(:,i)'];
    end
    scatter(subjDataX,subjDataY, 15, 'ko', 'filled')
    % Plot lines connecting
    for ii = 1:num_subj
        plot(1:size(curr_data,2),curr_data(ii,:), 'Color', subs_col((mod(ii,size(subs_col,1))+1),:))
    end
    if logRT
       yLimits=[6,7];
    else
       yLimits=[400,1050];
    end
end
ylim(yLimits);

hold off

%plot single items in each condition:
cond_shuffle=[7 9 8 10];
curr_data=curr_data(:,cond_shuffle);
num_subj=size(curr_data,1);
meanRT=nanmean(curr_data);
%compare sinalge items, v/nv:
subjAv=nanmean(curr_data,2);
withinEr=curr_data-repmat(subjAv,1,num_cond*2);
SEM=nanstd(withinEr)/sqrt(size(curr_data,1));

color=color(cond_shuffle,:);
    
subplot(1,2,2);
bar(1:length(meanRT),zeros(size(meanRT)),'FaceColor','none');
hold on
for i=1:num_cond*2%
    bar(i,meanRT(i),'FaceColor',color(i,:));
    errorbar(i,meanRT(i),SEM(i),'k');
end

if logRT
    ylabel('logRT','Fontsize',16);
    yLimits=[6.5,6.6];
else
    ylabel('RT(ms)','Fontsize',16);
    yLimits=[700,770];
end
title('simPEL11 v/nv block: violation','Fontsize',20);
set(gca,'XTickLabel',Xticks(cond_shuffle),'FontSize', 14,'XTickLabelRotation',25);

% Plot subject data
if plotSubj
    subjDataX=[];
    subjDataY=[];
    for i=1:num_cond
        subjDataX = [subjDataX ones(1,num_subj)*i];
        subjDataY = [subjDataY curr_data(:,i)'];
    end
    scatter(subjDataX,subjDataY, 15, 'ko', 'filled')
    % Plot lines connecting
    for ii = 1:num_subj
        plot(1:num_cond,curr_data(ii,1:num_cond), 'Color', subs_col((mod(ii,size(subs_col,1))+1),:))
    end
    if logRT
       yLimits=[6.3,7];
    else
       yLimits=[500,1050];
    end
end
ylim(yLimits);

hold off


col1=1;
col2=2;
fprintf('v/nv block: Brem: single items: v vs. nv: \n');
[h,p,ci,stats] = ttest(curr_data(:,col1),curr_data(:,col2));
p

col1=3;
col2=4;
fprintf('v/nv block: Bforg: single items: v vs. nv: \n');
[h,p,ci,stats] = ttest(curr_data(:,col1),curr_data(:,col2));
p

%remove one ss with nans:
curr_data=curr_data([1:17 19:end],:);
num_cond=size(curr_data,2)/2;
n=size(curr_data,1);
Y=reshape(curr_data,[numel(curr_data),1]);
S=repmat((1:n)',num_cond*2,1);
F1=[ones(n*num_cond,1);ones(n*num_cond,1)*2];%Brem/Bforg
F2=reshape(repmat(1:num_cond,n,2),size(Y));%violation/no-violation
fprintf('ANOVA single items: F1:AB memory: rem/forg, F2: violation/no-violation \n');
X=[Y,F1,F2,S];
RMAOV2_mod(X,alpha,showtable);

