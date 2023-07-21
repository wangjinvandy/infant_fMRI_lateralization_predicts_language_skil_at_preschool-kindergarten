%% First level analysis, written by Jin Wang 3/15/2019
% You should define your conditions, onsets, duration, TR.
% The repaired images will be deweighted from 1 to 0.01 in the first level
% estimation (we changed the art_redo.m, which uses art_deweight.txt as default to deweight the scans to art_redo_jin.txt, which we uses the art_repaired.txt to deweight scans).
% The difference between art_deiweghted.txt and art_repaired.txt is to the
% former one is more wide spread. It not only mark the scans which was repaired to be deweighted, but also scans around it to be deweighted.
% The 6 movement parameters we got from realignment is added into the model regressors to remove the small motion effects on data.

% Make sure you run clear all before running this code. This is to clear all existing data structure which might be left by previous analysis in the work space.

% This code is for ELP project specifically to deal with its repeated runs and run-<> is after acq- that would cause run-01 is after run-02 when specifying the model.

addpath(genpath('/n/gaab_mri_l3/Lab/DMC-Gaab2/data/Analysis_JinWang/fmri_code/')); %This is the code path
spm_path='/n/gaab_mri_l3/Lab/DMC-Gaab2/data/Analysis_JinWang/fmri_code/spm12/spm12_inf'; %This is your spm path
addpath(genpath(spm_path));

%define your data path
data=struct();
root='/n/gaab_mri_l3/Lab/DMC-Gaab2/data/Analysis_JinWang';  %your project path
subjects={};
data_info='/n/gaab_mri_l3/Lab/DMC-Gaab2/data/Analysis_JinWang/BabyBold_INF_code/sublist.xlsx';
if isempty(subjects)
    M=readtable(data_info);
    subjects=M.participant_all;
end

analysis_folder='analysis/analysis_INF'; % the name of your first level modeling folder
%model_deweight='deweight'; % the deweigthed modeling folder, it will be inside of your analysis folder
global CCN
CCN.preprocessed='BabyBold_sp'; % your data folder
CCN.session='INF'; % the time points you want to analyze
CCN.func_pattern='*SPEECHRHYME*'; % the name of your functional folders
CCN.file='swrtask*.nii'; % the name of your preprocessed data (4d)
CCN.rpfile='multiReg_*.mat'; %the movement files

%%define your task conditions, each run is a cell, be sure it follows the
%%sequence of the output from count_repaired_acc_rt.m. The task follows
%%alphabetic order based on its file name Gram is before Plaus task. 
conditions=[]; %start with an empty conditions. 
conditions{1}={'Forward' 'Backward' 'REST'};

%TR
TR=0.72; 

%define your contrasts, make sure your contrasts and your weights should be
%matched.
contrasts={'Forward_vs_Backward' 'Forward_vs_REST' 'Backward_vs_REST' 'Forward+Backward_vs_REST'};

%adjust the contrast by adding six 0s into the end of each session
weights={[1 -1 0]...
    [1 0 -1]...
    [0 1 -1]...
    [1 1 -2]};

%%%%%%%%%%%%%%%%%%%%%%%%Do not edit below here%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check if you define your contrasts in a correct way
if length(weights)~=length(contrasts)
    error('the contrasts and the weights are not matched');
end      

% Initialize
%addpath(spm_path);
spm('defaults','fmri');
spm_jobman('initcfg');
spm_figure('Create','Graphics','Graphics');

% Dependency and sanity checks
if verLessThan('matlab','R2013a')
    error('Matlab version is %s but R2013a or higher is required',version)
end

%Start to analyze the data from here
try
    for i=1:length(subjects)
        fprintf('work on subject %s', subjects{i});
        CCN.subject=[root '/' CCN.preprocessed '/' subjects{i}];
        %specify the outpath,create one if it does not exist
        out_path=[CCN.subject '/' analysis_folder];
        if ~exist(out_path)
            mkdir(out_path)
        end
        
        %find folders in func
        CCN.functional_dirs='[subject]/[session]/func/[func_pattern]/';
        functional_dirs=expand_path(CCN.functional_dirs);
        
                
        %load the functional data, 6 mv parameters, and event onsets
        mv=[];
        swfunc=[];
        P=[];
        onsets=[];
        dur=[];
        for j=1:length(functional_dirs)
            swfunc{j}=expand_path([functional_dirs{j} '[file]']);
%manually calculated onset times (in seconds) and duration base don the INF programs. 
            onsets{j}{1}=[0 37 91 127 235 271]; %forward
            dur{j}{1}=[18 18 18 18 18 18];
            onsets{j}{2}=[19 55 109 145 181 217]; %backward
            dur{j}{2}=[18 18 18 18 18 18];
            onsets{j}{3}=[73 163 199 253 289]; %rest1
            dur{j}{3}=[18 18 18 18 18];

            rp_file=expand_path([functional_dirs{j} '[rpfile]']);
            mv{j}=rp_file; 
        end

        data.swfunc=swfunc;
        
        
        %pass the experimental design information to data
        data.conditions=conditions;
        data.onsets=onsets;
        data.dur=dur;
        data.mv=mv;
        
        %run the firstlevel modeling and estimation (with deweighting)
        mat=firstlevel_4d_inf(data, out_path, TR);
  
        %run the contrasts
       
        contrast_f(mat,contrasts,weights);
        
    end
    
catch e
    rethrow(e)
    %display the errors
end