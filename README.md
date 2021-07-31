# simPEL_public

This includes the code for running the task and the analyses reported in "_Mnemonic prediction errors promote detailed memories_" https://psyarxiv.com/b5z8a/ Bein, Plotkin & Davachi Currently under revision.

Data and stimuli: https://osf.io/thzub/ 

simPEL9 = Experiment 1
simPEL11_associative = Experiment 2
simPEL11_item = Experiment 3

## running_exp_code folder:
This is the code we used to run the experiments.
Experiments 2 and 3 were run using exactly the same code, only the instructions to the participants were different.

##### These scripts are parallel in both Exp1 and Exp2/3, so I only explain re Exp1:
_simPEL9_randomize_images_: creates the sequences for each participant. We ran it before each participant was run.

_simPEL9_randomize_images_violateResp_: a version in which we tried to violate the response (in the ectual study we made sure not to violate the response) - we never run this study.

_simPEL9_learning_day1_: runs day1 of the learning task (predictive pairs learning)

_simPEL9_learning_day2_: runs day2 of the learning task (reminder, violation/no-violation phase).

_simPEL9_recognition_: runs the item recognition test.

_simPEL9_Explicit_: runs the associative memory test.

_do_math_distractor_: runs the math distractor (we did that between the item and the associative memory test).

The other files are auxilary functions, called by the other scripts.

* In the Exp2_3 folder there's also the script _simPEL11_randomize_images_for_sophie.m_. This one is identical to the other randomize_images code, but I commented this one a bit more for Sophie Nolden, from Yee Lee Shing's lab, who is runing a life span version of this.

## Analysis folder:

_simPEL9_analyse_AllLearningPhases_AB_mem_Oded_: analyse each participant's learning data (predictive pairs learning, reminder, violation phases), creates a matlab structure with the group data.

_simPEL9_group_learning_analyses_Oded_: run the group level analyses on the learning data.
In the paper - this is all reported in the Supplementary. The accuracy reported is based on this script. This script also includes some additional plotting and RT analysis, but the better way to analyse RTs is using gLMM (Lo 2015), so I did it this way, and this is reported in the paper, based on the R code (see below)

_simPEL9_analyse_recExplicit_Oded_:analyse each participant's item memory and associative memory data, creates a matlab structure with the group data.

_simPEL9_group_rec_explicit_analysesOded_: run the group lavel analyses for the item and associative memory data - analyses in this code are reported in the paper.

* Exp2 and 3 were the same only with the modification of the task. simPEL11associative... refers to Exp2, simPEL11item refers to Exp3.
* 
* rmANOVA folder: scripts for running rmANOVA in MATLAB, used in other analysis scripts.
* add here about the change accuracy.
I did the learning files, stoped at simPEL9_recExplicit..

