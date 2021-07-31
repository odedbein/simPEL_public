# simPEL_public

This includes the code for running the task and the analyses reported in "_Mnemonic prediction errors promote detailed memories_" https://psyarxiv.com/b5z8a/ Bein, Plotkin & Davachi Currently under revision.

Data and stimuli: https://osf.io/thzub/ 

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



