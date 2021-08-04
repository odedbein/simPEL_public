function simPEL9_11_compare_AB_conf(simPEL11_All_subs_explicit,simPEL9_All_subs_explicit)
%run;
%[simPEL11_All_subs_recognition, simPEL11_All_subs_explicit]=simPEL11associative_analyse_recExplicit_Oded()
%run:
%[simPEL9_All_subs_recognition, simPEL9_All_subs_explicit]=simPEL9_analyse_recExplicit_Oded();

%% gather simPEL9 (Experiment 1) data:
%good_subj based on low AB memory:
good_subj=[2 3 6 8:17 19:20 22 24:26 31:39]; %that is 28 subs
%choose_conds=[3 6 2 5]; 
Exp1_data=simPEL9_All_subs_explicit.Part2.CountEachRatingAvDist(good_subj,:);
%% gather simPEL11 (Experiment 2) data:
%good_subj based on low AB memory:
good_subj=[1:2 4:7 9:15 17:21 23:32];
Exp2_data=simPEL11_All_subs_explicit.Part2.CountEachRatingAvDist(good_subj,:);

%% plotting stuff:
color=[0        0.4470    0.7410
       0.8500   0.3250    0.0980
       0        0.4470    0.7410
       0.8500   0.3250    0.0980];

%% analyze count confidence rates:
color=[0        0.4470    0.7410
       0.8500   0.3250    0.0980
       0        0.4470    0.7410
       0.8500   0.3250    0.0980];
npairs=36;
data=Exp1_data;
%exclude confidence rates
Exc=1:5;
curr_data=[];
num_cond=2;
for i=1:num_cond*2
    curr_data=[curr_data sum(data(:,((i-1)*6+Exc(end)+1:i*6)),2)./npairs];
end
curr_data=curr_data(:,[1,3,2,4]);
conf_diff1=curr_data(:,[1,2])-curr_data(:,[3,4]);

data=Exp2_data;
curr_data=[];
num_cond=2;
for i=1:num_cond*2
    curr_data=[curr_data sum(data(:,((i-1)*6+Exc(end)+1:i*6)),2)./npairs];
end
curr_data=curr_data(:,[1,3,2,4]);
conf_diff2=curr_data(:,[1,2])-curr_data(:,[3,4]);


curr_data=[conf_diff2 conf_diff1];

meanAcc=mean(curr_data);
SEM=std(curr_data)/sqrt(size(curr_data,1));
f=figure;
set(f,'name','Exp1/2 conf_diff','numbertitle','off');
bar(1:length(meanAcc),meanAcc,'FaceColor',[1,1,1]);
hold on
for i=1:length(meanAcc) %
    bar(i,meanAcc(i),'FaceColor',color(i,:));
    errorbar(i,meanAcc(i),SEM(i),'k');
end
%ylim([0.8 4.5]);
Xticks={'EXP2: violation','EXP2: no-violation','EXP1: violation','EXP1: no-violation'};
set(gca,'XTickLabel',Xticks,'XTickLabelRotation',45)
set(gca, 'FontSize', 14)
title(sprintf('percentage confidence rates(hits-FA), Exclude 1-%d',Exc(end)),'FontSize', 18)

hold off

disp('rates confidence');
fprintf('v: exp1/2 t-tests \n');
[h,p,ci,stats]=ttest(curr_data(:,1),curr_data(:,3));
p

fprintf('nv: exp1/2 t-tests \n');
[h,p,ci,stats]=ttest(curr_data(:,2),curr_data(:,4));
p

n=size(curr_data,1);
data=reshape(curr_data,numel(curr_data),1);
data=[data ([ones(n*2,1);ones(n*2,1)*2])];
data=[data ([ones(n,1);ones(n,1)*2;ones(n,1);ones(n,1)*2])];
data=[data ([1:n 1:n 1:n 1:n])'];
fprintf('Repeated measures ANOVA on average confidence rates, V1 is EXPERIMENT(Exp2/Exp1), V2 is condition (v/nv)\n\n');
RMAOV2_mod(data);


%% compute differences per confidence levels:
npairs=36;
Exp1_diff=[Exp1_data(:,1:6)-Exp1_data(:,7:12) Exp1_data(:,13:18)-Exp1_data(:,19:24)]/npairs;
Exp1_diff=reshape(Exp1_diff,[28,6,2]);
Exp1_diff=mean(Exp1_diff,3);
mean1_diff=mean(Exp1_diff);
Exp2_diff=[Exp2_data(:,1:6)-Exp2_data(:,7:12) Exp2_data(:,13:18)-Exp2_data(:,19:24)]/npairs;
Exp2_diff=reshape(Exp2_diff,[28,6,2]);
Exp2_diff=mean(Exp2_diff,3);
mean2_diff=mean(Exp2_diff);

x = [1 2 3 4 5 6];
figure;
hold on
bar(x,mean1_diff,.8,'FaceColor',[0 0 0])
bar(x,mean2_diff,.5,'FaceColor',[0.7 0.7 0.7]) %'None','EdgeColor',
legend({'Exp1','Exp2'},'Location','northwest','FontSize',14)
xlim([.5 6.5]);
title('percentage confidence rates(hits-misses)','FontSize', 18)
xlabel('confidence level','FontSize', 14)
hold off

fprintf('Conf1 exp1/2 t-tests \n');
[h,p,ci,stats]=ttest2(Exp1_diff(:,1),Exp2_diff(:,1));
p

fprintf('conf6 exp1/2 t-tests \n');
[h,p,ci,stats]=ttest2(Exp1_diff(:,6),Exp2_diff(:,6));
p

fprintf('conf6-conf1 exp1/2 t-tests \n');
[h,p,ci,stats]=ttest2(Exp1_diff(:,6)-Exp1_diff(:,1),Exp2_diff(:,6)-Exp2_diff(:,1));
p



