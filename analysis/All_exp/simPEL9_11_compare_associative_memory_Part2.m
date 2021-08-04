%this script is not meant to run as a function - but by executing sections
engram=0;
if engram
    mydir='/data/Bein';
else
    mydir='/Volumes/data/Bein';
end

project_dir=fullfile(mydir,'/simPEL/simPEL11_onlyConsSameSimilarItemTask');
analysis_dir=fullfile(project_dir,'analysis');
results_dir=fullfile(analysis_dir,'files_forR');
if ~exist(results_dir); mkdir(results_dir); end



%% compare associative memory for the original predictive pair (Part 2), between study1 and study3, and btw study 2 and study 3:
% 1. run the scripts to create simPEL9_All_subs_explicit
good_subj1=[2 3 6 8:17 19:20 22 24:26 31:39]; %that is 28 subs as reported in my paper
Exp1_data=simPEL9_All_subs_explicit.Part2.OldNewAvDist(good_subj1,[1,3]);
Exp1_col=reshape(Exp1_data,numel(Exp1_data),1);

% 2. run the scripts to create simPEL11_All_subs_explicit (associative)
good_subj2=[1:2 4:7 9:15 17:21 23:32];%that is 28 subs as reported in my paper
Exp2_data=simPEL11_All_subs_explicit.Part2.OldNewAvDist(good_subj2,[1,3]);
Exp2_col=reshape(Exp2_data,numel(Exp2_data),1);

% 3. run the scripts to create simPEL11_All_subs_explicit (item)
%here we take all participants
good_subj3=1:28;
Exp3_data=simPEL11_All_subs_explicit.Part2.OldNewAvDist(good_subj3,[1,3]);
Exp3_col=reshape(Exp3_data,numel(Exp3_data),1);

%make the model:
% experiment participant Bmem: rem/forg viol: violation/No-violation ret_trial_type:old/sim/new response: old/sim/new
%some of these work just bc we had identical N in both studies, change if not. 
n=length(good_subj1);
subj=[repmat(good_subj1',2,1);repmat((good_subj2'+200),2,1);repmat((good_subj3'+300),2,1)];%response,violation/no-violation, rem/forg,sim/old item
Experiment=[ones(length(Exp1_col),1);ones(length(Exp2_col),1)*2;ones(length(Exp3_col),1)*3];
viol_intact=[ones(length(Exp1_col)/2,1);zeros(length(Exp1_col)/2,1);ones(length(Exp2_col)/2,1);zeros(length(Exp2_col)/2,1);ones(length(Exp3_col)/2,1);zeros(length(Exp3_col)/2,1)];

data=[Exp1_col;Exp2_col;Exp3_col];
model_data=[Experiment subj viol_intact data];

header={'Experiment' 'subject' 'viol_intact' 'Bmem_rates'};
% set the table:
T=array2table(model_data);
T.Properties.VariableNames=header;
%write it up:
results_fname='simPEL9_11ass_11item_ExpPart2.txt';
filename=fullfile(results_dir,results_fname);
writetable(T,filename,'Delimiter','\t')