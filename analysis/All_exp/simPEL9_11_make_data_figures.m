function simPEL9_11_make_data_figures(simPEL11_All_subs_recognition,simPEL9_All_subs_recognition)
%run;
%[simPEL11_All_subs_recognition, simPEL11_All_subs_explicit]=simPEL11associative_analyse_recExplicit_Oded()
%run:
%[simPEL9_All_subs_recognition, simPEL9_All_subs_explicit]=simPEL9_analyse_recExplicit_Oded();

%% gather simPEL9 (Experiment 1) data:
%good_subj based on low AB memory:
good_subj=[2 3 6 8:17 19:20 22 24:26 31:39]; %that is 28 subs
choose_conds=[3 6 2 5]; %'Brem: violation-old: old' 'Brem: no-violation-old: old','Brem: violation-old: similar','Brem: no-violation-old: similar'
Exp1_data=simPEL9_All_subs_recognition.All.Bmem.response_typePerCond(good_subj,choose_conds); %remove "new" responses, change the order of similar/old

%% gather simPEL11 (Experiment 2) data:
%good_subj based on low AB memory:
good_subj=[1:2 4:7 9:15 17:21 23:32];
Exp2_data=simPEL11_All_subs_recognition.All.Bmem.response_typePerCond(good_subj,choose_conds); %remove "new" responses, change the order of similar/old


%% gather simPEL11 item (Experiment 3) data (used it to make the data figure for the talk in Yee Lee Shing's lab Oct. 2020:
subjects_numbers=[1:15 18:30];

good_subj=1:numel(subjects_numbers); %[1:9 11:14 16];
Exp3_data=simPEL11_All_subs_recognition.All.Bmem.response_typePerCond(good_subj,choose_conds); %remove "new" responses, change the order of similar/old


%% plot:

color=[
    0.4500    0.0250    0.0980
    0         0.1470    0.3410
    0.6500    0.1250    0.0980
    0         0.2470    0.5410
    ];

color=[
    0.6500    0.1250    0.0980
    0         0.2470    0.5410
    0.8500    0.3250    0.0980
    0         0.4470    0.7410
    ];

color=[
    0.6500    0.1250    0.0980
    0         0.2470    0.5410
    0.6500    0.1250    0.0980
    0         0.2470    0.5410
    ];
figure;
for curr_exp=1:2
    if curr_exp==1
        curr_data=Exp1_data;
        exp_title='Experiment 1';
        
    elseif curr_exp==2
        curr_data=Exp2_data;
        exp_title='Experiment 2';
    end
    
    
    num_cond=size(curr_data,2);
    meanAcc=nanmean(curr_data);
    SEM=nanstd(curr_data)./sqrt(size(curr_data,1));
    
    subplot(1,2,curr_exp);
    hold on
    bar(1:num_cond,zeros(1,num_cond),'FaceColor','none');
    for i=1:num_cond %
        bar(i,meanAcc(i),'FaceColor',color(i,:));
        errorbar(i,meanAcc(i),SEM(i),'k');
    end
    ylim([0 0.9]);
    xlim([0.2 num_cond+0.8]);
    ylabel('Response rates','Fontsize',16);
    title(exp_title,'Fontsize',24);

    hold off
    
end

%% make the AB forgotten figure for the supplementary
%good_subj based on low AB memory:
good_subj=[2 3 6 8:17 19:20 22 24:26 31:39]; %that is 28 subs
choose_conds=[15 18 14 17 21 24 20 23];
Exp1_data=simPEL9_All_subs_recognition.All.Bmem.response_typePerCond(good_subj,choose_conds); %remove "new" responses, change the order of similar/old

%% gather simPEL11 (Experiment 2) data:
%good_subj based on low AB memory:
good_subj=[1:2 4:7 9:15 17:21 23:32];
Exp2_data=simPEL11_All_subs_recognition.All.Bmem.response_typePerCond(good_subj,choose_conds); %remove "new" responses, change the order of similar/old


color=[
    0.6500    0.1250    0.0980
    0         0.2470    0.5410
    0.8500    0.3250    0.0980
    0         0.4470    0.7410
    
    0.6500    0.1250    0.0980
    0         0.2470    0.5410
    0.8500    0.3250    0.0980
    0         0.4470    0.7410
    ];

figure;
for curr_exp=1:2
    if curr_exp==1
        curr_data=Exp1_data;
        exp_title='Experiment 1';
        
    elseif curr_exp==2
        curr_data=Exp2_data;
        exp_title='Experiment 2';
    end
    
    
    num_cond=size(curr_data,2);
    meanAcc=nanmean(curr_data);
    SEM=nanstd(curr_data)./sqrt(size(curr_data,1));
    
    subplot(2,1,curr_exp);
    hold on
    bar(1:num_cond,zeros(1,num_cond),'FaceColor','none');
    for i=1:num_cond %
        bar(i,meanAcc(i),'FaceColor',color(i,:));
        errorbar(i,meanAcc(i),SEM(i),'k');
    end
    ylim([0 0.9]);
    xlim([0.2 num_cond+0.8]);
    ylabel('Response rates','Fontsize',16);
    title(exp_title,'Fontsize',24);

    hold off
    
end




end