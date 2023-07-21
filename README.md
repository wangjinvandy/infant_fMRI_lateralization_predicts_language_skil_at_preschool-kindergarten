# infant_fMRI_lateralization_predicts_language_skil_at_preschool-kindergarten
This is the code used for the paper "Left-lateralization of the superior temporal gyrus during speech processing in sleeping infants predicts language skills in kindergarten: a task-based fMRI study. "   
#^ The preprocessing of infant fMRI data can be found https://github.com/TeddyTuresky/infant_fmri_preproc, written by Ted Turesky.    
#^ firstlevel_generate_inf.m is the firstlevel code. It relies on spm12 (we renamed it as spm12_inf in the code) toolbox, where we modified the spm_hrf.m to make it an infant hrf.  
#^ secondlevel_generate.m is the secondlevel code to generate group activation map using one-sample t test. This code also relies on spm12 toolbox.   
#^ top_voxel_mask.sh is the code used to extract the top activated voxels. This code relies on AFNI.  
#^ getbetas.m is the code used to extract the brain activation from each individualized top activated voxels. This code needs the marsbar toolbox.   
#^ inf_fmri_data_longitudinal.xlsx is a file where all the brain measures, behavioral measures, and participants used for this current study.   
