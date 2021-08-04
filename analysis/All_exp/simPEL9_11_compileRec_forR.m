function simPEL9_11_compileRec_forR(simPEL11_All_subs_recognition,simPEL9_All_subs_recognition,engram)
%run;
%[simPEL11_All_subs_recognition, simPEL11_All_subs_explicit,simPEL11_All_subs_explicitDependency]=simPEL11associative_analyse_recExplicit_Oded()
%run:
%[simPEL9_All_subs_recognition, simPEL9_All_subs_explicit,simPEL9_All_subs_explicitDependency]=simPEL9_analyse_recExplicit_Oded();
% use the reported_analyses version.

if engram
    mydir='/data/Bein';
else
    mydir='/Volumes/data/Bein';
end

project_dir=fullfile(mydir,'/simPEL/simPEL11_onlyConsSameSimilarItemTask');
analysis_dir=fullfile(project_dir,'analysis');
results_dir=fullfile(analysis_dir,'files_forR');
if ~exist(results_dir); mkdir(results_dir); end


%choose and organize rec conditions:
%not like the data, but like I'll read it - that mataches the choose_cond
%columns
%that's useful to build the model - but then I don't use it.
header.rec={'Brem: violation-old: old','Brem: violation-old: similar','Brem: violation-old: new', 'Brem: no-violation-old: old','Brem: no-violation-old: similar','Brem: no-violation-old: new',...
            'Brem: violation-sim: old','Brem: violation-sim: similar','Brem: violation-sim: new', 'Brem: no-violation-sim: old','Brem: no-violation-sim: similar','Brem: no-violation-sim: new'...
            'Bforg: violation-old: old','Bforg: violation-old: similar','Bforg: violation-old: new', 'Bforg: no-violation-old: old','Bforg: no-violation-old: similar','Bforg: no-violation-old: new',...
            'Bforg: violation-sim: old','Bforg: violation-sim: similar', 'Bforg: violation-sim: new', 'Bforg: no-violation-sim: old','Bforg: no-violation-sim: similar','Bforg: no-violation-sim: new'...
             };
choose_conds=[3 2 1 6 5 4 9 8 7 12 11 10 15 14 13 18 17 16 21 20 19 24 23 22];
%% gather simPEL9 (Experiment 1) data:
%good_subj based on low AB memory:
good_subj1=[2 3 6 8:17 19:20 22 24:26 31:39]; %that is 28 subs as reported in my paper

Exp1_data=simPEL9_All_subs_recognition.All.Bmem.response_typePerCond(good_subj1,choose_conds); %remove "new" responses, change the order of similar/old
Exp1_col=reshape(Exp1_data,numel(Exp1_data),1);

%% gather simPEL11 (Experiment 2) data:
%good_subj based on low AB memory:
good_subj2=[1:2 4:7 9:15 17:21 23:32];%that is 28 subs as reported in my paper
Exp2_data=simPEL11_All_subs_recognition.All.Bmem.response_typePerCond(good_subj2,choose_conds); %remove "new" responses, change the order of similar/old
Exp2_col=reshape(Exp2_data,numel(Exp2_data),1);

%make the model:
% experiment participant Bmem: rem/forg viol: violation/No-violation ret_trial_type:old/sim/new response: old/sim/new
%some of these work just bc we had identical N in both studies, change if not. 
n=length(good_subj1);
subj=[repmat(good_subj1',2*2*2*3,1);repmat((good_subj2'+200),2*2*2*3,1)];%response,violation/no-violation, rem/forg,sim/old item
Experiment=[ones(length(Exp1_col),1);ones(length(Exp1_col),1)*2];
Bmem=[ones(length(Exp1_col)/2,1);zeros(length(Exp1_col)/2,1);ones(length(Exp2_col)/2,1);zeros(length(Exp2_col)/2,1)];
ret_trial_type=repmat([ones(n*6,1);zeros(n*6,1)],2*2,1); %rem/forg,Experiments
viol=repmat([ones(n*3,1);zeros(n*3,1)],2*2*2,1); %rem/forg,sim/old item,Experiments
response=repmat([ones(n,1)*2;ones(n,1);zeros(n,1)],2*2*2*2,1); %violation/no-violation, rem/forg,sim/old item,Experiments

data=[Exp1_col;Exp2_col];
model_data=[Experiment subj Bmem ret_trial_type viol response data];

header={'Experiment' 'subject' 'Bmemory' 'ret_trial_type' 'viol' 'response' 'rates'};
% set the table:
T=array2table(model_data);
T.Properties.VariableNames=header;
%write it up:
results_fname='simPEL9_11_recognition.txt';
filename=fullfile(results_dir,results_fname);
writetable(T,filename,'Delimiter','\t')

end