%% this script is not meant to run as a function - but by executing sections
%set stats parameters:
showtable=1;
alpha=0.05;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% recognition - analyse accuracy:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% the begining is w/o looking at AB remembered/forgotten separately, so some of it is not reported in the paper.
%% The results figure in the main paper is based on the script in the All_exp folder

plotSubj=0;
header.rec={'violation old','no-violation old','violation sim-lure','no-violation sim-lure','new'};
subjects_numbers=[2 4 5:8 12:18 20:30 34 36:39 41 43:48 80:82]';

%bad memory Exp Part 2 - less than 40%:
%overall bad rates in both conditions below 40%:
%subj num: 2,6,7,12,24,27,29,37,38,39,41,
%subj num 1 is excluded due to technical error.

%good_subj based on low AB memory:
good_subj=[2 3 6 8:17 19:20 22 24:26 31:39]; %that is 28 subs

%analyze accuracy:
data=simPEL9_All_subs_recognition.All.accuracy_rates;
header.rec(1,2:end+1)=header.rec(1,1:end);
header.rec{1,1}='subjects';
AccStruct.rec=header.rec;
AccStruct.rec(2:length(subjects_numbers)+1,1)=num2cell(subjects_numbers);
AccStruct.rec(2:length(subjects_numbers)+1,2:end)=num2cell(data);

data=simPEL9_All_subs_recognition.All.accuracy_rates(good_subj,:);
curr_data=data(:,:);%use if want to remove

num_subj=size(curr_data,1);
num_cond=size(curr_data,2);
meanAcc=mean(curr_data);
subjAv=mean(curr_data,2);
withinEr=curr_data-repmat(subjAv,1,num_cond);
SEM=std(withinEr)/sqrt(size(curr_data,1));

color=[
       0         0.2470    0.5410
       0.6500    0.1250    0.0980
       
       0         0.4470    0.7410
       0.8500    0.3250    0.0980
       
       0.7       0.7       0.7
       
       ];

color=[color;color]; %for both halves
figure;
subs_col=get(gca,'colororder');

hold on
bar(1:num_cond,meanAcc,'FaceColor',[1,1,1]);
for i=1:length(meanAcc) %
    bar(i,meanAcc(i),'FaceColor',color(i,:));
    errorbar(i,meanAcc(i),SEM(i),'k');
end
ylim([0 1.2]);
xlim([0.2 num_cond+0.8]);
Xticks=header.rec(2:end);
set(gca,'XTickLabel',Xticks,'XTickLabelRotation',25)
set(gca, 'FontSize', 14)
ylabel('accuracy','Fontsize',16);   
title('simPEL9 Recognition','Fontsize',24);

% Plot subject data
if plotSubj
    subjDataX=[];
    subjDataY=[];
    for i=1:num_cond
        subjDataX = [subjDataX ones(1,num_subj)*i];
        subjDataY = [subjDataY curr_data(:,i)'];
    end
    scatter(subjDataX,subjDataY, 15, 'ko', 'filled')
    % Plot lines connecting within/across
    for ii = 1:num_subj
        plot(1:num_cond,curr_data(ii,1:num_cond), 'Color', subs_col((mod(ii,size(subs_col,1))+1),:))
    end
    
end
hold off

%% recognition - analyse all response types
plotSubj=0;

%good_subj based on low AB memory (same as the sample above):
good_subj=[2 3 6 8:17 19:20 22 24:26 31:39]; %that is 28 subs

%not like the data, but like I'll read it - changing the similar/old
%columns
header.rec={'violation-old: old','violation-old: similar', 'no-violation-old: old','no-violation-old: similar',...
            'violation-sim: old','violation-sim: similar', 'no-violation-sim: old','no-violation-sim: similar',...
            'new: old','new: similar'...
            };


%analyze accuracy:
data=simPEL9_All_subs_recognition.All.response_typePerCond(:,[3 2 6 5 9 8 12 11 15 14]); %remove "new" responses, change the order of similar/old
header.rec(1,2:end+1)=header.rec(1,1:end);
header.rec{1,1}='subjects';
AccStruct.rec=header.rec;
AccStruct.rec(2:length(subjects_numbers)+1,1)=num2cell(subjects_numbers');
AccStruct.rec(2:length(subjects_numbers)+1,2:end)=num2cell(data);

data=data(good_subj,:);
%data=data([1:24 26:end],:);%use if want to remove participant with perfect AB memory
curr_data=data(:,:);%use if want to remove conditions - did it before

num_subj=size(curr_data,1);
num_cond=size(curr_data,2);
meanAcc=mean(curr_data);
subjAv=mean(curr_data,2);
withinEr=curr_data-repmat(subjAv,1,num_cond);
SEM=std(withinEr)/sqrt(size(curr_data,1));

color=[
       0         0.1470    0.3410
       0         0.2470    0.5410
       
       0.4500    0.0250    0.0980
       0.6500    0.1250    0.0980
       
       0         0.4470    0.7410
       0         0.6470    0.9410
       
       0.8500    0.3250    0.0980
       0.9500    0.5250    0.2980
       
       0.7       0.7       0.7
       0.9       0.9       0.9
       
       ];

color=[color;color]; %for both halves
figure;
subplot(2,1,1);
subs_col=get(gca,'colororder');
hold on
bar(1:num_cond,meanAcc,'FaceColor',[1,1,1]);
for i=1:length(meanAcc) %
    bar(i,meanAcc(i),'FaceColor',color(i,:));
    errorbar(i,meanAcc(i),SEM(i),'k');
end
ylim([0 1]);
xlim([0.2 num_cond+0.8]);
Xticks=header.rec(2:end);
set(gca,'XTick', 1:num_cond, 'XTickLabel',Xticks,'XTickLabelRotation',25)
set(gca, 'FontSize', 14)
ylabel('accuracy','Fontsize',16);   
title('simPEL9 by halves - all response types','Fontsize',24);

% Plot subject data
if plotSubj
    subjDataX=[];
    subjDataY=[];
    for i=1:num_cond
        subjDataX = [subjDataX ones(1,num_subj)*i];
        subjDataY = [subjDataY curr_data(:,i)'];
    end
    scatter(subjDataX,subjDataY, 15, 'ko', 'filled')
    % Plot lines connecting within/across
    for ii = 1:num_subj
        plot(1:num_cond,curr_data(ii,1:num_cond), 'Color', subs_col((mod(ii,size(subs_col,1))+1),:))
    end
    
end
hold off

% now plot it by grouping old/similaer
shuffle_conds=[1 3 2 4 5 7 6 8 9 10];
curr_data=data(:,shuffle_conds);%use if want to remove conditions - did it before
num_subj=size(curr_data,1);
num_cond=size(curr_data,2);
meanAcc=mean(curr_data);
subjAv=mean(curr_data,2);
withinEr=curr_data-repmat(subjAv,1,num_cond);
SEM=std(withinEr)/sqrt(size(curr_data,1));

color=color(shuffle_conds,:);
figure;
%subplot(2,1,2);
subs_col=get(gca,'colororder');
hold on
bar(1:num_cond,meanAcc,'FaceColor',[1,1,1]);
for i=1:length(meanAcc) %
    bar(i,meanAcc(i),'FaceColor',color(i,:));
    errorbar(i,meanAcc(i),SEM(i),'k');
end
ylim([0 1]);
xlim([0.2 num_cond+0.8]);
Xticks=Xticks(shuffle_conds);
set(gca,'XTick', 1:num_cond, 'XTickLabel',Xticks,'XTickLabelRotation',25)
set(gca, 'FontSize', 14)
ylabel('accuracy','Fontsize',16);   
title('simPEL9 - similar/old responses','Fontsize',24);

% Plot subject data
if plotSubj
    subjDataX=[];
    subjDataY=[];
    for i=1:num_cond
        subjDataX = [subjDataX ones(1,num_subj)*i];
        subjDataY = [subjDataY curr_data(:,i)'];
    end
    scatter(subjDataX,subjDataY, 15, 'ko', 'filled')
    % Plot lines connecting within/across
    for ii = 1:num_subj
        for cc=1:2:num_cond
        plot(cc:cc+1,curr_data(ii,cc:cc+1), 'Color', subs_col((mod(ii,size(subs_col,1))+1),:))
        end
    end
    
end
hold off

%% stats - average on violation and no-violation, just to know whether subjects differntiate
%% identical old from similar lures:
data=[mean(curr_data(:,1:2),2) mean(curr_data(:,3:4),2) mean(curr_data(:,5:6),2) mean(curr_data(:,7:8),2) curr_data(:,9:10)]; %save for later use
stats_data=data(:,1:4);

fprintf('old items, old responses: M = %.2f, SD: %.2f\n',mean(stats_data(:,1)),std(stats_data(:,1)));
fprintf('old items, similar responses: M = %.2f, SD: %.2f\n',mean(stats_data(:,2)),std(stats_data(:,2)));
fprintf('similar items, old responses: M = %.2f, SD: %.2f\n',mean(stats_data(:,3)),std(stats_data(:,3)));
fprintf('similar items, similar responses: M = %.2f, SD: %.2f\n',mean(stats_data(:,4)),std(stats_data(:,4)));
fprintf('new items, old responses: M = %.2f, SD: %.2f\n',mean(data(:,5)),std(data(:,5)));
fprintf('new items, similar responses: M = %.2f, SD: %.2f\n',mean(data(:,6)),std(data(:,6)));


fprintf('old items: old vs. similar responses: \n');
[h,p,ci,stats] = ttest(stats_data(:,1),stats_data(:,2))
fprintf('similar items: old vs. similar responses: \n');
[h,p,ci,stats] = ttest(stats_data(:,3),stats_data(:,4))
fprintf('old items: old vs. foils: old \n');
[h,p,ci,stats] = ttest(data(:,1),data(:,5))
fprintf('similar items: similar vs. foils: similar \n');
[h,p,ci,stats] = ttest(data(:,4),data(:,6))



num_cond=size(stats_data,2)/2;
n=size(stats_data,1);
Y=reshape(stats_data,size(stats_data,1)*size(stats_data,2),1);
S=repmat([1:n]',num_cond*2,1);
F1=[ones(n*num_cond,1);ones(n*num_cond,1)*2];%old/similar items
F2=reshape(repmat(1:num_cond,n,2),size(Y));%old/similar responses
fprintf('ANOVA old items: F1:item statuse: old/similar, F2: old/similar responses \n');
X=[Y,F1,F2,S];
RMAOV2_mod(X,alpha,showtable);

%% stats:
data=curr_data; %save for later use
%within old items - old vs. similar responses
curr_data=data(:,1:4);
num_cond=size(curr_data,2)/2;
n=size(curr_data,1);
Y=reshape(curr_data,size(curr_data,1)*size(curr_data,2),1);
S=repmat([1:n]',num_cond*2,1);
F1=[ones(n*num_cond,1);ones(n*num_cond,1)*2];%old/similar
F2=reshape(repmat(1:num_cond,n,2),size(Y));%violation/no-violation
display(sprintf('ANOVA old items: F1:response type: old/similar, F2: violation/no-violation \n'));
X=[Y,F1,F2,S];
RMAOV2_mod(X,alpha,showtable);


%t-tests for simple effects:
col1=1;
col2=2;
display(sprintf('old: v/nv \n'));
[h,p,ci,stats] = ttest(curr_data(:,col1),curr_data(:,col2));
p

col1=3;
col2=4;
display(sprintf('similar: v/nv  \n'));
[h,p,ci,stats] = ttest(curr_data(:,col1),curr_data(:,col2));
p

%within old responses for old vs. similar items
curr_data=data(:,[1:2 5:6]);
num_cond=size(curr_data,2)/2;
n=size(curr_data,1);
Y=reshape(curr_data,size(curr_data,1)*size(curr_data,2),1);
S=repmat([1:n]',num_cond*2,1);
F1=[ones(n*num_cond,1);ones(n*num_cond,1)*2];%old/similar
F2=reshape(repmat(1:num_cond,n,2),size(Y));%violation/no-violation
display(sprintf('ANOVA old responses: F1:item type: old/similar, F2: violation/no-violation \n'));
X=[Y,F1,F2,S];
RMAOV2_mod(X,alpha,showtable);

%within similar responses for old vs. similar items
stats_data=data(:,[3:4 7:8]);
num_cond=size(stats_data,2)/2;
n=size(stats_data,1);
Y=reshape(stats_data,size(stats_data,1)*size(stats_data,2),1);
S=repmat([1:n]',num_cond*2,1);
F1=[ones(n*num_cond,1);ones(n*num_cond,1)*2];%old/similar
F2=reshape(repmat(1:num_cond,n,2),size(Y));%violation/no-violation
fprintf('ANOVA similar responses: F1:item type: old/similar, F2: violation/no-violation \n');
X=[Y,F1,F2,S];
RMAOV2_mod(X,alpha,showtable);


fprintf('similar lures, similar responses: v/nv \n');
col1=3;
col2=4;
[~,p,ci,stats] = ttest(stats_data(:,col1),stats_data(:,col2));
CohenD=nanmean(stats_data(:,col1)-stats_data(:,col2))/nanstd(stats_data(:,col1)-stats_data(:,col2));
fprintf('t: %.2f, p: %.3f, d: %.2f \n',stats.tstat,p,CohenD);


%within similar items, old vs. similar responses
stats_data=data(:,[5:6 7:8]);
num_cond=size(stats_data,2)/2;
n=size(stats_data,1);
Y=reshape(stats_data,size(stats_data,1)*size(stats_data,2),1);
S=repmat([1:n]',num_cond*2,1);
F1=[ones(n*num_cond,1);ones(n*num_cond,1)*2];%old/similar
F2=reshape(repmat(1:num_cond,n,2),size(Y));%violation/no-violation
fprintf('ANOVA similar items: F1:responses: old/similar, F2: violation/no-violation \n');
X=[Y,F1,F2,S];
RMAOV2_mod(X,alpha,showtable);


fprintf('similar items, old responses: v/nv \n');
col1=1;
col2=2;
[~,p,ci,stats] = ttest(stats_data(:,col1),stats_data(:,col2));
CohenD=nanmean(stats_data(:,col1)-stats_data(:,col2))/nanstd(stats_data(:,col1)-stats_data(:,col2));
fprintf('t: %.2f, p: %.3f, d: %.2f \n',stats.tstat,p,CohenD);


%% recognition - analyse all response types separated by B item memory - remembered vs. forgotten B items:
plotSubj=0;

%not like the data, but like I'll read it - changing the similar/old
%columns
header.rec={'Brem: violation-old: old','Brem: violation-old: similar', 'Brem: no-violation-old: old','Brem: no-violation-old: similar',...
            'Brem: violation-sim: old','Brem: violation-sim: similar', 'Brem: no-violation-sim: old','Brem: no-violation-sim: similar',...
            'Bforg: violation-old: old','Bforg: violation-old: similar', 'Bforg: no-violation-old: old','Bforg: no-violation-old: similar',...
            'Bforg: violation-sim: old','Bforg: violation-sim: similar', 'Bforg: no-violation-sim: old','Bforg: no-violation-sim: similar',...
             };

%analyze accuracy:
choose_conds=[3 2 6 5 9 8 12 11 15 14 18 17 21 20 24 23];
data=simPEL9_All_subs_recognition.All.Bmem.response_typePerCond(:,choose_conds); %remove "new" responses, change the order of similar/old
data=data(good_subj,:);

%data=data([1:24 26:end],:); %removing this participant here if want to exclude AB perfect memory
header.rec(1,2:end+1)=header.rec(1,1:end);
header.rec{1,1}='subjects';
Xticks=header.rec(2:end);

color=[
       0         0.1470    0.3410
       0         0.2470    0.5410
       
       0.4500    0.0250    0.0980
       0.6500    0.1250    0.0980
       
       0         0.4470    0.7410
       0         0.6470    0.9410
       
       0.8500    0.3250    0.0980
       0.9500    0.5250    0.2980
       
       ];

color=repmat(color,2,1); %for rem/forg,


%plot it by grouping old/similaer
shuffle_conds=[1 3 9 11 2 4 10 12 5 7 13 15 6 8 14 16];
curr_data=data(:,shuffle_conds);%use if want to remove conditions - did it before
num_subj=size(curr_data,1);
num_cond=size(curr_data,2);
meanAcc=nanmean(curr_data);
subjAv=nanmean(curr_data,2);
withinEr=curr_data-repmat(subjAv,1,num_cond);
SEM=nanstd(withinEr)/sqrt(size(curr_data,1));
color=color(shuffle_conds,:);


figure;
subs_col=get(gca,'colororder');
hold on
bar(1:num_cond,meanAcc(1:num_cond),'FaceColor',[1,1,1]);
for i=1:num_cond %
    bar(i,meanAcc(i),'FaceColor',color(i,:));
    errorbar(i,meanAcc(i),SEM(i),'k');
end
ylim([0 1]);
xlim([0.2 num_cond+0.8]);
Xticks=Xticks(shuffle_conds);
set(gca,'XTick', 1:num_cond, 'XTickLabel',Xticks,'XTickLabelRotation',25)
set(gca, 'FontSize', 14)
ylabel('accuracy','Fontsize',16);   
title('simPEL9 - similar/old responses','Fontsize',24);

% Plot subject data
if plotSubj
    subjDataX=[];
    subjDataY=[];
    for i=1:num_cond
        subjDataX = [subjDataX ones(1,num_subj)*i];
        subjDataY = [subjDataY curr_data(:,i)'];
    end
    scatter(subjDataX,subjDataY, 15, 'ko', 'filled')
    % Plot lines connecting within/across
    for ii = 1:num_subj
        for cc=1:2:num_cond
        plot(cc:cc+1,curr_data(ii,cc:cc+1), 'Color', subs_col((mod(ii,size(subs_col,1))+1),:))
        end
    end
    
end
hold off

%% stats:
%old responses vs. similar responses on old items, only rem AB associations
stats_data=curr_data(:,[1:2 5:6]);

fprintf('Brem: V: old items, old responses: M = %.2f, SD: %.2f\n',mean(stats_data(:,1)),std(stats_data(:,1)));
fprintf('Brem: NV: old items, old responses: M = %.2f, SD: %.2f\n',mean(stats_data(:,2)),std(stats_data(:,2)));
fprintf('Brem: V: old items, similar responses: M = %.2f, SD: %.2f\n',mean(stats_data(:,3)),std(stats_data(:,3)));
fprintf('Brem: NV: old items, similar responses: M = %.2f, SD: %.2f\n \n',mean(stats_data(:,4)),std(stats_data(:,4)));


fprintf('Brem: old items: old responses: v/nv \n');
col1=1;
col2=2;
[~,p,ci,stats] = ttest(stats_data(:,col1),stats_data(:,col2));
CohenD=nanmean(stats_data(:,col1)-stats_data(:,col2))/nanstd(stats_data(:,col1)-stats_data(:,col2));
fprintf('t: %.2f, p: %.3f, d: %.2f \n',stats.tstat,p,CohenD);

fprintf('Brem: old items: similar responses: v/nv: \n');
col1=3;
col2=4;
[~,p,ci,stats] = ttest(stats_data(:,col1),stats_data(:,col2));
CohenD=nanmean(stats_data(:,col1)-stats_data(:,col2))/nanstd(stats_data(:,col1)-stats_data(:,col2));
fprintf('t: %.2f, p: %.3f, d: %.2f \n',stats.tstat,p,CohenD);

%ANOVA: B-rem, old items: F1:response type: old/similar, F2: violation/no-violation
num_cond=size(stats_data,2)/2;
n=size(stats_data,1);
Y=reshape(stats_data,size(stats_data,1)*size(stats_data,2),1);
S=repmat([1:n]',num_cond*2,1);
F1=[ones(n*num_cond,1);ones(n*num_cond,1)*2];%old/similar
F2=reshape(repmat(1:num_cond,n,2),size(Y));%violation/no-violation
fprintf('ANOVA B-rem, old items: F1:response type: old/similar, F2: violation/no-violation \n');
X=[Y,F1,F2,S];
RMAOV2_mod_oded(X,alpha,showtable);

%old responses to old items vs. similar items, only rem AB associations
stats_data=curr_data(:,[1:2 9:10]);
fprintf('Brem: V: similar items, old responses: M = %.2f, SD: %.2f\n',mean(stats_data(:,3)),std(stats_data(:,3)));
fprintf('Brem: NV: similar items, old responses: M = %.2f, SD: %.2f\n \n',mean(stats_data(:,4)),std(stats_data(:,4)));
col1=3;
col2=4;
[~,p,ci,stats] = ttest(stats_data(:,col1),stats_data(:,col2));
CohenD=nanmean(stats_data(:,col1)-stats_data(:,col2))/nanstd(stats_data(:,col1)-stats_data(:,col2));
fprintf('t: %.2f, p: %.3f, d: %.2f \n',stats.tstat,p,CohenD);

%ANOVA: B-rem, old responses. F1: item: old/similar, F2: violation/no-violation
num_cond=size(stats_data,2)/2;
n=size(stats_data,1);
Y=reshape(stats_data,size(stats_data,1)*size(stats_data,2),1);
S=repmat([1:n]',num_cond*2,1);
F1=[ones(n*num_cond,1);ones(n*num_cond,1)*2];%old/similar
F2=reshape(repmat(1:num_cond,n,2),size(Y));%violation/no-violation
fprintf('ANOVA B-rem, old responses: F1:item: old/similar, F2: violation/no-violation \n');
X=[Y,F1,F2,S];
RMAOV2_mod_oded(X,alpha,showtable);

%similar responses to old items vs. similar items, only rem AB associations
stats_data=curr_data(:,[5:6 13:14]);
fprintf('Brem: V: old items, similar responses: M = %.2f, SD: %.2f\n',mean(stats_data(:,1)),std(stats_data(:,1)));
fprintf('Brem: NV: old items, similar responses: M = %.2f, SD: %.2f\n \n',mean(stats_data(:,2)),std(stats_data(:,2)));

fprintf('Brem: V: similar items, similar responses: M = %.2f, SD: %.2f\n',mean(stats_data(:,3)),std(stats_data(:,3)));
fprintf('Brem: NV: similar items, similar responses: M = %.2f, SD: %.2f\n \n',mean(stats_data(:,4)),std(stats_data(:,4)));

fprintf('Brem: similar items: similar responses: v/nv \n');
col1=3;
col2=4;
[~,p,ci,stats] = ttest(stats_data(:,col1),stats_data(:,col2));
CohenD=nanmean(stats_data(:,col1)-stats_data(:,col2))/nanstd(stats_data(:,col1)-stats_data(:,col2));
fprintf('t: %.2f, p: %.3f, d: %.2f \n',stats.tstat,p,CohenD);


fprintf('Brem: old items: similar responses: v/nv \n');
col1=1;
col2=2;
[~,p,ci,stats] = ttest(stats_data(:,col1),stats_data(:,col2));
CohenD=nanmean(stats_data(:,col1)-stats_data(:,col2))/nanstd(stats_data(:,col1)-stats_data(:,col2));
fprintf('t: %.2f, p: %.3f, d: %.2f \n',stats.tstat,p,CohenD);

%ANOVA: B-rem, similar responses. F1: item: old/similar, F2: violation/no-violation
num_cond=size(stats_data,2)/2;
n=size(stats_data,1);
Y=reshape(stats_data,size(stats_data,1)*size(stats_data,2),1);
S=repmat([1:n]',num_cond*2,1);
F1=[ones(n*num_cond,1);ones(n*num_cond,1)*2];%old/similar
F2=reshape(repmat(1:num_cond,n,2),size(Y));%violation/no-violation
fprintf('ANOVA B-rem, similar responses: F1:item: old/similar item, F2: violation/no-violation \n');
X=[Y,F1,F2,S];
RMAOV2_mod_oded(X,alpha,showtable);


%FORG AB: old responses vs. similar responses on old items
% subj num 48: 100% mem for b items - remove for anova with forgotten
%old responses on old items, by rem/forg AB associations
stats_data=curr_data([1:24 26:end],[3:4 7:8]);

fprintf('Bforg: V: old items, old responses: M = %.2f, SD: %.2f\n',mean(stats_data(:,1)),std(stats_data(:,1)));
fprintf('Bforg: NV: old items, old responses: M = %.2f, SD: %.2f\n',mean(stats_data(:,2)),std(stats_data(:,2)));
fprintf('Bforg: V: old items, similar responses: M = %.2f, SD: %.2f\n',mean(stats_data(:,3)),std(stats_data(:,3)));
fprintf('Bforg: NV: old items, similar responses: M = %.2f, SD: %.2f\n \n',mean(stats_data(:,4)),std(stats_data(:,4)));


fprintf('Bforg: old items: old responses: v/nv \n');
col1=1;
col2=2;
[~,p,ci,stats] = ttest(stats_data(:,col1),stats_data(:,col2));
CohenD=nanmean(stats_data(:,col1)-stats_data(:,col2))/nanstd(stats_data(:,col1)-stats_data(:,col2));
fprintf('t: %.2f, p: %.3f, d: %.2f \n',stats.tstat,p,CohenD);

fprintf('Bforg: old items: similar responses: v/nv: \n');
col1=3;
col2=4;
[~,p,ci,stats] = ttest(stats_data(:,col1),stats_data(:,col2));
CohenD=nanmean(stats_data(:,col1)-stats_data(:,col2))/nanstd(stats_data(:,col1)-stats_data(:,col2));
fprintf('t: %.2f, p: %.3f, d: %.2f \n',stats.tstat,p,CohenD);

%ANOVA: B-forg, old items: F1:response type: old/similar, F2: violation/no-violation
num_cond=size(stats_data,2)/2;
n=size(stats_data,1);
Y=reshape(stats_data,size(stats_data,1)*size(stats_data,2),1);
S=repmat([1:n]',num_cond*2,1);
F1=[ones(n*num_cond,1);ones(n*num_cond,1)*2];%old/similar
F2=reshape(repmat(1:num_cond,n,2),size(Y));%violation/no-violation
fprintf('ANOVA B-rem, old items: F1:response type: old/similar, F2: violation/no-violation \n');
X=[Y,F1,F2,S];
RMAOV2_mod_oded(X,alpha,showtable);

%for violation identical old items, Brem/Bforg, old/similar responses:
% subj num 48: 100% mem for b items - remove for anova with forgotten
stats_data=curr_data([1:24 26:end],[1,3,5,7]);
fprintf('Bforg-Brem: old items: old responses: \n');
col1=1;
col2=2;
[~,p,ci,stats] = ttest(stats_data(:,col1),stats_data(:,col2));
CohenD=nanmean(stats_data(:,col1)-stats_data(:,col2))/nanstd(stats_data(:,col1)-stats_data(:,col2));
fprintf('t: %.2f, p: %.3f, d: %.2f \n',stats.tstat,p,CohenD);

fprintf('Bforg-Brem: old items: similar responses: \n');
col1=3;
col2=4;
[~,p,ci,stats] = ttest(stats_data(:,col1),stats_data(:,col2));
CohenD=nanmean(stats_data(:,col1)-stats_data(:,col2))/nanstd(stats_data(:,col1)-stats_data(:,col2));
fprintf('t: %.2f, p: %.3f, d: %.2f \n',stats.tstat,p,CohenD);
%ANOVA: old responses on old items, by rem/forg AB associations
num_cond=size(stats_data,2)/2;
n=size(stats_data,1);
Y=reshape(stats_data,size(stats_data,1)*size(stats_data,2),1);
S=repmat([1:n]',num_cond*2,1);
F1=[ones(n*num_cond,1);ones(n*num_cond,1)*2];%'old'/'similar'
F2=reshape(repmat(1:num_cond,n,2),size(Y));%B-rem/B-forg
fprintf('ANOVA violation old items: F1:response type:old/similar, F2: B-rem/B-forg \n');
X=[Y,F1,F2,S];
RMAOV2_mod_oded(X,alpha,showtable);

%% recognition - analyse all response types separated by B item memory don't plot forgotten items:
plotSubj=0;
%good_subj based on low AB memory (same as the sample above):
good_subj=[2 3 6 8:17 19:20 22 24:26 31:39]; %that is 28 subs

%not like the data, but like I'll read it - changing the similar/old
%columns
header.rec={'Brem: violation-old: old','Brem: violation-old: similar', 'Brem: no-violation-old: old','Brem: no-violation-old: similar',...
            'Brem: violation-sim: old','Brem: violation-sim: similar', 'Brem: no-violation-sim: old','Brem: no-violation-sim: similar',...
             };

%analyze accuracy:
choose_conds=[3 2 6 5 9 8 12 11];
data=simPEL9_All_subs_recognition.All.Bmem.response_typePerCond(:,choose_conds); %remove "new" responses, change the order of similar/old
data=data(good_subj,:);
header.rec(1,2:end+1)=header.rec(1,1:end);
header.rec{1,1}='subjects';
Xticks=header.rec(2:end);

color=[
       0         0.1470    0.3410
       0         0.2470    0.5410
       
       0.4500    0.0250    0.0980
       0.6500    0.1250    0.0980
       
       0         0.4470    0.7410
       0         0.6470    0.9410
       
       0.8500    0.3250    0.0980
       0.9500    0.5250    0.2980
       
       ];

%plot it by grouping old/similaer
shuffle_conds=[1 3 2 4 5 7 6 8];
curr_data=data(:,shuffle_conds);%use if want to remove conditions - did it before
num_subj=size(curr_data,1);
num_cond=size(curr_data,2);
meanAcc=nanmean(curr_data);
subjAv=nanmean(curr_data,2);
withinEr=curr_data-repmat(subjAv,1,num_cond);
SEM=nanstd(withinEr)/sqrt(size(curr_data,1));
color=color(shuffle_conds,:);


figure;
%subplot(2,1,2);
subs_col=get(gca,'colororder');
hold on
bar(1:num_cond,meanAcc(1:num_cond),'FaceColor',[1,1,1]);
for i=1:num_cond %
    bar(i,meanAcc(i),'FaceColor',color(i,:));
    errorbar(i,meanAcc(i),SEM(i),'k');
end
ylim([0 1]);
xlim([0.2 num_cond+0.8]);
Xticks=Xticks(shuffle_conds);
set(gca,'XTick', 1:num_cond, 'XTickLabel',Xticks,'XTickLabelRotation',25)
set(gca, 'FontSize', 14)
ylabel('accuracy','Fontsize',16);   
title('simPEL9 - similar/old responses','Fontsize',24);

% Plot subject data
if plotSubj
    subjDataX=[];
    subjDataY=[];
    for i=1:num_cond
        subjDataX = [subjDataX ones(1,num_subj)*i];
        subjDataY = [subjDataY curr_data(:,i)'];
    end
    scatter(subjDataX,subjDataY, 15, 'ko', 'filled')
    % Plot lines connecting within/across
    for ii = 1:num_subj
        for cc=1:2:num_cond
        plot(cc:cc+1,curr_data(ii,cc:cc+1), 'Color', subs_col((mod(ii,size(subs_col,1))+1),:))
        end
    end
    
end
hold off

%t-tests for simple effects:
col1=1;
col2=2;
fprintf('B-rem, old items: old responses: v/nv \n');
[h,p,ci,stats] = ttest(curr_data(:,col1),curr_data(:,col2));
p

col1=3;
col2=4;
fprintf('B-forg: similar: v/nv  \n');
[h,p,ci,stats] = ttest(curr_data(:,col1),curr_data(:,col2));
p


%old responses vs. similar responses on old items, only rem AB associations
stats_data=curr_data(:,1:4);
num_cond=size(stats_data,2)/2;
n=size(stats_data,1);
Y=reshape(stats_data,size(stats_data,1)*size(stats_data,2),1);
S=repmat([1:n]',num_cond*2,1);
F1=[ones(n*num_cond,1);ones(n*num_cond,1)*2];%old/similar
F2=reshape(repmat(1:num_cond,n,2),size(Y));%violation/no-violation
fprintf('ANOVA old items: F1:response type: old/similar, F2: violation/no-violation \n');
X=[Y,F1,F2,S];
RMAOV2_mod(X,alpha,showtable);

% calculate and plot the "discrimination score" - saying old on old vs. saying similar on old.
disc_score=[curr_data(:,1)-curr_data(:,3) curr_data(:,2)-curr_data(:,4)];
%ttest:
col1=1;
col2=2;
fprintf('B-rem, old items, discrimination score (saying "old" - saying "similar"): v/nv \n');
[h,p,ci,stats] = ttest(disc_score(:,col1),disc_score(:,col2));
p

num_subj=size(disc_score,1);
num_cond=size(disc_score,2);
meanAcc=nanmean(disc_score);
subjAv=nanmean(disc_score,2);
withinEr=disc_score-repmat(subjAv,1,num_cond);
SEM=nanstd(withinEr)/sqrt(size(disc_score,1));
color=color([1 2],:);


figure;
%subplot(2,1,2);
subs_col=get(gca,'colororder');
hold on
bar(1:num_cond,meanAcc(1:num_cond),'FaceColor',[1,1,1]);
for i=1:num_cond %
    bar(i,meanAcc(i),'FaceColor',color(i,:));
    errorbar(i,meanAcc(i),SEM(i),'k');
end
ylim([0.2 0.6]);
xlim([0.2 num_cond+0.8]);
Xticks={'violation','no-violation'};
set(gca,'XTick', 1:num_cond, 'XTickLabel',Xticks,'XTickLabelRotation',25)
set(gca, 'FontSize', 14)
ylabel('discremination score (old items, responded "old"-"similar")','Fontsize',16);   
title('only B-remembered','Fontsize',24);
hold off

%% recognition - INCLUDE NEW RESPONSES analyse all response types separated by B item memory - remembered vs. forgotten B items:
plotSubj=0;

%not like the data, but like I'll read it - changing the similar/old
%columns
header.rec={'Brem: violation-old: old','Brem: violation-old: similar','Brem: violation-old: new', 'Brem: no-violation-old: old','Brem: no-violation-old: similar','Brem: no-violation-old: new',...
            'Brem: violation-sim: old','Brem: violation-sim: similar','Brem: violation-sim: new', 'Brem: no-violation-sim: old','Brem: no-violation-sim: similar','Brem: no-violation-sim: new'...
            'Bforg: violation-old: old','Bforg: violation-old: similar','Bforg: violation-old: new', 'Bforg: no-violation-old: old','Bforg: no-violation-old: similar','Bforg: no-violation-old: new',...
            'Bforg: violation-sim: old','Bforg: violation-sim: similar', 'Bforg: violation-sim: new', 'Bforg: no-violation-sim: old','Bforg: no-violation-sim: similar','Bforg: no-violation-sim: new'...
             };

%analyze accuracy:
choose_conds=[3 2 1 6 5 4 9 8 7 12 11 10 15 14 13 18 17 16 21 20 19 24 23 22];
data=simPEL9_All_subs_recognition.All.Bmem.response_typePerCond(:,choose_conds); %remove "new" responses, change the order of similar/old
header.rec(1,2:end+1)=header.rec(1,1:end);
header.rec{1,1}='subjects';
AccStruct.rec=header.rec;
AccStruct.rec(2:length(subjects_numbers)+1,1)=num2cell(subjects_numbers');
AccStruct.rec(2:length(subjects_numbers)+1,2:end)=num2cell(data);

data=data(good_subj,:);
%data=data([1:17 19:end],:);
curr_data=data(:,:);%use if want to remove conditions - did it before
num_subj=size(curr_data,1);
num_cond=size(curr_data,2);
meanAcc=nanmean(curr_data);
% subjAv=nanmean(curr_data,2);
% withinEr=curr_data-repmat(subjAv,1,num_cond);
SEM=nanstd(curr_data)/sqrt(size(curr_data,1));

color=[
       0         0.1470    0.3410
       0         0.2470    0.5410
       1         1         1
       0.4500    0.0250    0.0980
       0.6500    0.1250    0.0980
       1         1         1
       0         0.4470    0.7410
       0         0.6470    0.9410
       1         1         1
       0.8500    0.3250    0.0980
       0.9500    0.5250    0.2980
       1         1         1
       ];

color=repmat(color,2,1); %for rem/forg,

figure;
%subplot(2,1,1);
subs_col=get(gca,'colororder');
hold on
bar(1:num_cond,meanAcc(1:num_cond),'FaceColor',[1,1,1]);
for i=1:num_cond %
    bar(i,meanAcc(i),'FaceColor',color(i,:));
    errorbar(i,meanAcc(i),SEM(i),'k');
end
ylim([0 1]);
xlim([0.2 num_cond+0.8]);
Xticks=header.rec(2:end);
set(gca,'XTick', 1:num_cond, 'XTickLabel',Xticks,'XTickLabelRotation',25)
set(gca, 'FontSize', 14)
ylabel('accuracy','Fontsize',16);   
title('simPEL9 - all response types','Fontsize',24);

% Plot subject data
if plotSubj
    subjDataX=[];
    subjDataY=[];
    for i=1:num_cond
        subjDataX = [subjDataX ones(1,num_subj)*i];
        subjDataY = [subjDataY curr_data(:,i)'];
    end
    scatter(subjDataX,subjDataY, 15, 'ko', 'filled')
    % Plot lines connecting within/across
    for ii = 1:num_subj
        plot(1:num_cond,curr_data(ii,1:num_cond), 'Color', subs_col((mod(ii,size(subs_col,1))+1),:))
    end
    
end
hold off

%now plot it by grouping old/similaer
shuffle_conds=[1 4 2 5 3 6 7 10 8 11 9 12 [1 4 2 5 3 6 7 10 8 11 9 12]+12];
curr_data=data(:,shuffle_conds);%use if want to remove conditions - did it before
num_subj=size(curr_data,1);
num_cond=size(curr_data,2);
meanAcc=nanmean(curr_data);
SEM=nanstd(curr_data)/sqrt(size(curr_data,1));

color=color(shuffle_conds,:);


figure;
subs_col=get(gca,'colororder');
hold on
bar(1:num_cond,meanAcc(1:num_cond),'FaceColor',[1,1,1]);
for i=1:num_cond %
    bar(i,meanAcc(i),'FaceColor',color(i,:));
    errorbar(i,meanAcc(i),SEM(i),'k');
end
ylim([0 1]);
xlim([0.2 num_cond+0.8]);
Xticks=Xticks(shuffle_conds);
set(gca,'XTick', 1:num_cond, 'XTickLabel',Xticks,'XTickLabelRotation',25)
set(gca, 'FontSize', 14)
ylabel('accuracy','Fontsize',16);   
title('simPEL9 - similar/old responses','Fontsize',24);

% Plot subject data
if plotSubj
    subjDataX=[];
    subjDataY=[];
    for i=1:num_cond
        subjDataX = [subjDataX ones(1,num_subj)*i];
        subjDataY = [subjDataY curr_data(:,i)'];
    end
    scatter(subjDataX,subjDataY, 15, 'ko', 'filled')
    % Plot lines connecting within/across
    for ii = 1:num_subj
        for cc=1:2:num_cond
        plot(cc:cc+1,curr_data(ii,cc:cc+1), 'Color', subs_col((mod(ii,size(subs_col,1))+1),:))
        end
    end
    
end
hold off

% plot only Brem:
curr_conds=1:12;

figure;
%subplot(2,1,2);
subs_col=get(gca,'colororder');
hold on
bar(curr_conds,meanAcc(curr_conds),'FaceColor',[1,1,1]);
for i=curr_conds %
    bar(i,meanAcc(i),'FaceColor',color(i,:));
    errorbar(i,meanAcc(i),SEM(i),'k');
end
ylim([0 1]);
xlim([0.2 curr_conds(end)+0.8]);
Xticks=Xticks(curr_conds);
set(gca,'XTick', curr_conds, 'XTickLabel',Xticks,'XTickLabelRotation',25)
set(gca, 'FontSize', 14)
ylabel('accuracy','Fontsize',16);   
title('simPEL9 - similar/old responses','Fontsize',24);

% Plot subject data
if plotSubj
    subjDataX=[];
    subjDataY=[];
    for i=curr_conds
        subjDataX = [subjDataX ones(1,num_subj)*i];
        subjDataY = [subjDataY curr_data(:,i)'];
    end
    scatter(subjDataX,subjDataY, 15, 'ko', 'filled')
    % Plot lines connecting within/across
    for ii = 1:num_subj
        for cc=1:2:curr_conds(end)
        plot(cc:cc+1,curr_data(ii,cc:cc+1), 'Color', subs_col((mod(ii,size(subs_col,1))+1),:))
        end
    end
    
end
hold off

% "new responses":
stats_data=curr_data(:,[5,6,11,12]);
fprintf('Brem: V: old items, new responses: M = %.2f, SD: %.2f\n',mean(stats_data(:,1)),std(stats_data(:,1)));
fprintf('Brem: NV: old items, new responses: M = %.2f, SD: %.2f\n \n',mean(stats_data(:,2)),std(stats_data(:,2)));

fprintf('Brem: old items: new responses: v/nv \n');
col1=1;
col2=2;
[~,p,ci,stats] = ttest(stats_data(:,col1),stats_data(:,col2));
CohenD=nanmean(stats_data(:,col1)-stats_data(:,col2))/nanstd(stats_data(:,col1)-stats_data(:,col2));
fprintf('t: %.2f, p: %.3f, d: %.2f \n',stats.tstat,p,CohenD);

% "new responses": 
fprintf('Brem: V: similar items, new responses: M = %.2f, SD: %.2f\n',mean(stats_data(:,3)),std(stats_data(:,3)));
fprintf('Brem: NV: similar items, new responses: M = %.2f, SD: %.2f\n \n',mean(stats_data(:,4)),std(stats_data(:,4)));

fprintf('Brem: similar items: new responses: v/nv \n');
col1=3;
col2=4;
[~,p,ci,stats] = ttest(stats_data(:,col1),stats_data(:,col2));
CohenD=nanmean(stats_data(:,col1)-stats_data(:,col2))/nanstd(stats_data(:,col1)-stats_data(:,col2));
fprintf('t: %.2f, p: %.3f, d: %.2f \n',stats.tstat,p,CohenD);


% "similar responses to similar items":
stats_data=curr_data(:,[9,10]);
fprintf('Brem: V: similar items, similar responses: M = %.2f, SD: %.2f\n',mean(stats_data(:,1)),std(stats_data(:,1)));
fprintf('Brem: NV: similar items, similar responses: M = %.2f, SD: %.2f\n \n',mean(stats_data(:,2)),std(stats_data(:,2)));

fprintf('Brem: similar items: similar responses: v/nv \n');
col1=1;
col2=2;
[~,p,ci,stats] = ttest(stats_data(:,col1),stats_data(:,col2));
CohenD=nanmean(stats_data(:,col1)-stats_data(:,col2))/nanstd(stats_data(:,col1)-stats_data(:,col2));
fprintf('t: %.2f, p: %.3f, d: %.2f \n',stats.tstat,p,CohenD);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% analyze explicit - part1 - AB memory
plotSubj=0;
%good_subj based on low AB memory (same as the sample above):
good_subj=[2 3 6 8:17 19:20 22 24:26 31:39]; %that is 28 subs

%analyze old/new
data=simPEL9_All_subs_explicit.Part1.OldNewAvDist;
header.Exp.part1={'violation','violation foil','no-violation','no-violation foil'};
%meanExplicit=mean(data);
header.Exp.part1(1,2:end+1)=header.Exp.part1(1,1:end);
header.Exp.part1{1,1}='subjects';
ExpStruct.part1=header.Exp.part1;
ExpStruct.part1(2:length(subjects_numbers)+1,1)=num2cell(subjects_numbers);
ExpStruct.part1(2:length(subjects_numbers)+1,2:end)=num2cell(data);

data=simPEL9_All_subs_explicit.Part1.OldNewAvDist(good_subj,:);
color=[0    0.4470    0.7410
    0.8500    0.3250    0.0980];

curr_data=data(:,[1,3]);
num_subj=size(curr_data,1);
num_cond=size(curr_data,2);
meanAcc=mean(curr_data);
subjAv=mean(curr_data,2);
withinEr=curr_data-repmat(subjAv,1,num_cond);
SEM=std(withinEr)/sqrt(size(curr_data,1));
f=figure;
set(f,'name','simPEL9','numbertitle','off');
bar(1:num_cond,meanAcc,'FaceColor',[1,1,1]);
hold on
for i=1:length(meanAcc) %
    bar(i,meanAcc(i),'FaceColor',color(i,:));
    errorbar(i,meanAcc(i),SEM(i),'k');
end
ylim([0.2 0.5]);
Xticks={'violation','no-violation'};
set(gca,'XTickLabel',Xticks)
set(gca, 'FontSize', 14)
title('% hits','FontSize', 18)

% Plot subject data
if plotSubj
    subjDataX = [ones(1,num_subj), ones(1,num_subj)*2];
    scatter(subjDataX,[curr_data(:,1);curr_data(:,2)]', 15, 'ko', 'filled')
    % Plot lines connecting rPlus/rMinus and assocPlus/assocMinus
    for ii = 1:num_subj
        plot([1,2], [curr_data(ii,1), curr_data(ii,2)])
    end
end
hold off

display(sprintf('t-tests AX/BX memory \n'));
display(sprintf('v/nv t-tests \n'));
[h,p,ci,stats]=ttest(curr_data(:,1),curr_data(:,2))

chance=ones(size(curr_data,1),1)*0.33;
display(sprintf('v diff from chance \n'));
[h,p,ci,stats]=ttest(curr_data(:,1),chance)

display(sprintf('nv diff from chance \n'));
[h,p,ci,stats]=ttest(curr_data(:,2),chance)


%% analyze average confidence rates:
plotSubj=0;
data=simPEL9_All_subs_explicit.Part1.AvConfRatesAvDist(good_subj,:);
color=[0        0.4470    0.7410
       0.8500   0.3250    0.0980
       0.8      0.8       0.8
       0.8      0.8       0.8];

curr_data=data(:,[1,3,2,4]);
num_subj=size(curr_data,1);
num_cond=size(curr_data,2);
meanAcc=mean(curr_data);
subjAv=mean(curr_data,2);
withinEr=curr_data-repmat(subjAv,1,num_cond);
SEM=std(withinEr)/sqrt(size(curr_data,1));
f=figure;
set(f,'name','simPEL9','numbertitle','off');
bar(1:length(meanAcc),meanAcc,'FaceColor',[1,1,1]);
hold on
for i=1:length(meanAcc) %
    bar(i,meanAcc(i),'FaceColor',color(i,:));
    errorbar(i,meanAcc(i),SEM(i),'k');
end
ylim([2.5 3.5]);
Xticks={'violation-hits','no-violation-hits','violation-foils','no-violation-foils'};
set(gca,'XTickLabel',Xticks,'XTickLabelRotation',25)
set(gca, 'FontSize', 14)
title('average confidence','FontSize', 18)

% Plot subject data
if plotSubj
subjDataX=[];
subjDataY=[];
for i=1:length(meanAcc)
    subjDataX = [subjDataX ones(1,num_subj)*i];
    subjDataY = [subjDataY curr_data(:,i)'];
end
scatter(subjDataX,subjDataY, 15, 'ko', 'filled')

% Plot lines connecting violation and no-violation
for ii = 1:num_subj
   plot([1,2], [curr_data(ii,1), curr_data(ii,2)]) 
   plot([3,4], [curr_data(ii,3), curr_data(ii,4)])
end
end
hold off

disp('average confidence rates');
fprintf('v/nv t-tests \n');
[h,p,ci,stats]=ttest(curr_data(:,1),curr_data(:,2));
p
fprintf('v/nv foils t-tests \n');
[h,p,ci,stats]=ttest(curr_data(:,3),curr_data(:,4));
p

fprintf('v-foils \n');
[h,p,ci,stats]=ttest(curr_data(:,1),curr_data(:,3));
p
fprintf('nv-foils \n');
[h,p,ci,stats]=ttest(curr_data(:,2),curr_data(:,4));
p

n=size(curr_data,1);
data=reshape(curr_data,numel(curr_data),1);
data=[data ([ones(n*2,1);ones(n*2,1)*2])];
data=[data ([ones(n,1);ones(n,1)*2;ones(n,1);ones(n,1)*2])];
data=[data ([1:n 1:n 1:n 1:n])'];
fprintf('Repeated measures ANOVA on average confidence rates, V1 is target(target/foils), V2 is condition (v/nv)\n\n');
RMAOV2_mod(data);


%% analyze count confidence rates:
data=simPEL9_All_subs_explicit.Part1.CountEachRatingAvDist(good_subj,:);
color=[0        0.4470    0.7410
       0.8500   0.3250    0.0980
       0.8      0.8       0.8
       0.8      0.8       0.8];
%exclude confidence rates
Exc=1:2;
curr_data=[];
num_cond=2;
for i=1:num_cond*2
    curr_data=[curr_data sum(data(:,((i-1)*6+Exc(end)+1:i*6)),2)];
end
curr_data=curr_data(:,[1,3,2,4]);
meanAcc=mean(curr_data);
subjAv=mean(curr_data,2);
withinEr=curr_data-repmat(subjAv,1,num_cond*2);
SEM=std(withinEr)/sqrt(size(curr_data,1));
f=figure;
set(f,'name','simPEL9','numbertitle','off');
bar(1:length(meanAcc),meanAcc,'FaceColor',[1,1,1]);
hold on
for i=1:length(meanAcc) %
    bar(i,meanAcc(i),'FaceColor',color(i,:));
    errorbar(i,meanAcc(i),SEM(i),'k');
end
%ylim([0.8 4.5]);
Xticks={'violation-hits','no-violation-hits','violation-foils','no-violation-foils'};
set(gca,'XTickLabel',Xticks,'XTickLabelRotation',45)
set(gca, 'FontSize', 14)
title(sprintf('num confidence rates, Exclude 1-%d',Exc(end)),'FontSize', 18)

% Plot subject data
subjDataX=[];
subjDataY=[];
for i=1:length(meanAcc)
    subjDataX = [subjDataX ones(1,num_subj)*i];
    subjDataY = [subjDataY curr_data(:,i)'];
end
scatter(subjDataX,subjDataY, 15, 'ko', 'filled')

% Plot lines connecting violation and no-violation
for ii = 1:num_subj
   plot([1,2], [curr_data(ii,1), curr_data(ii,2)]) 
   plot([3,4], [curr_data(ii,3), curr_data(ii,4)])
end
hold off

disp('count confidence');
fprintf('v/nv t-tests \n');
[h,p,ci,stats]=ttest(curr_data(:,1),curr_data(:,2));
p
fprintf('v/nv foils t-tests \n');
[h,p,ci,stats]=ttest(curr_data(:,3),curr_data(:,4));
p

fprintf('v-foils \n');
[h,p,ci,stats]=ttest(curr_data(:,1),curr_data(:,3));
p
fprintf('nv-foils \n');
[h,p,ci,stats]=ttest(curr_data(:,2),curr_data(:,4));
p

n=size(curr_data,1);
data=reshape(curr_data,numel(curr_data),1);
data=[data ([ones(n*2,1);ones(n*2,1)*2])];
data=[data ([ones(n,1);ones(n,1)*2;ones(n,1);ones(n,1)*2])];
data=[data ([1:n 1:n 1:n 1:n])'];
fprintf('Repeated measures ANOVA on average confidence rates, V1 is target(target/foils), V2 is condition (v/nv)\n\n');
RMAOV2_mod(data);

%% analyze explicit - PART2 - memory for the novel associations
%good_subj based on low AB memory (same as the sample above):
good_subj=[2 3 6 8:17 19:20 22 24:26 31:39]; %that is 28 subs

plotSubj=0;
%analyze old/new
data=simPEL9_All_subs_explicit.Part2.OldNewAvDist(good_subj,:);
header.Exp.part2={'AB-violation','violation foil','AB-no-violation','no-violation foil'};
%meanExplicit=mean(data);
header.Exp.part2(1,2:end+1)=header.Exp.part2(1,1:end);
header.Exp.part2{1,1}='subjects';
ExpStruct.part2=header.Exp.part2;
ExpStruct.part2(2:length(good_subj)+1,1)=num2cell(subjects_numbers(good_subj)');
ExpStruct.part2(2:length(good_subj)+1,2:end)=num2cell(data);

color=[0    0.4470    0.7410
    0.8500    0.3250    0.0980];

curr_data=data(:,[1,3]);
num_subj=size(curr_data,1);
num_cond=size(curr_data,2);
meanAcc=mean(curr_data);
SD_curr_data=std(curr_data);
subjAv=mean(curr_data,2);
withinEr=curr_data-repmat(subjAv,1,num_cond);
SEM=std(withinEr)/sqrt(size(curr_data,1));
f=figure;
set(f,'name','simPEL9','numbertitle','off');
bar(1:num_cond,meanAcc,'FaceColor',[1,1,1]);
hold on
for i=1:length(meanAcc) %
    bar(i,meanAcc(i),'FaceColor',color(i,:));
    errorbar(i,meanAcc(i),SEM(i),'k');
end
ylim([0.33 0.7]);
Xticks={'AB-violation','AB-no-violation'};
set(gca,'XTickLabel',Xticks)
set(gca, 'FontSize', 14)
title('AB memory: % hits','FontSize', 18)

% Plot subject data
if plotSubj
    subjDataX = [ones(1,num_subj), ones(1,num_subj)*2];
    scatter(subjDataX,[curr_data(:,1);curr_data(:,2)]', 15, 'ko', 'filled')
    % Plot lines connecting rPlus/rMinus and assocPlus/assocMinus
    for ii = 1:num_subj
        plot([1,2], [curr_data(ii,1), curr_data(ii,2)])
    end
end
hold off

fprintf('violation: M: %.2f, SD: %.2f \n', meanAcc(1),SD_curr_data(1));
fprintf('no-violation: M: %.2f, SD: %.2f \n', meanAcc(2),SD_curr_data(2));
fprintf('v/nv t-tests \n');
[h,p,ci,stats]=ttest(curr_data(:,1),curr_data(:,2))


