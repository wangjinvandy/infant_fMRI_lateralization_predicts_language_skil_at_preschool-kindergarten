%% Second level analysis
%currently I only write the one sample t test and
%regression analysis, which are the two mostly used.
%Jin Wang 3/19/2019


addpath(genpath('/n/gaab_mri_l3/Lab/DMC-Gaab2/data/Analysis_JinWang/fmri_code/')); %This is the code path
spm_path='/n/gaab_mri_l3/Lab/DMC-Gaab2/data/Analysis_JinWang/fmri_code/spm12/spm12_inf'; %This is your spm path
addpath(genpath(spm_path));

%define your data path
data=struct();
root='/n/gaab_mri_l3/Lab/DMC-Gaab2/data/Analysis_JinWang/BabyBold_sp';  %your project path
subjects={};
data_info='/n/gaab_mri_l3/Lab/DMC-Gaab2/data/Analysis_JinWang/BabyBold_INF_code/sublist.txt';
if isempty(subjects)
    M=readtable(data_info);
    subjects=M.participant_all;
end

out_dir='/n/gaab_mri_l3/Lab/DMC-Gaab2/data/Analysis_JinWang/BabyBold_sp/sp_secondlevel_INF_longitudinal'; % your output folder of secondlevel analysis results

% find the data. If you follow the data structure and naming convention in
% the walkthrough, you won't need to change this part
global CCN;
CCN.session='INF';
CCN.func_pattern='swrtask*';
analysis_folder='analysis/analysis_INF';

% choose your analysis method
test=1; %1 one-sample t test, 2 mutiple regression analysis, 3 two-sample t test

% if you have covariates in your second level modeling
cov=0; % 1 if you have covariates or 0 if you do not have covariates

% Ignore these lines if you do not have covariates.
% You do not need to comment the following information out if you put cov=0
% because it won't be read in the analysis. Only when you put
%cov=1, the following covariates will be read in the analysis.

% define your covariates of control for your one-sample t test. Or define your covariates of interest for your multiple regression
% analysis if you have covariates.
cov_num=2; %number of your covariates
%you can define as many covariates as you want by adding name and values.
if cov==1
    name=[];
    name{1}=''; % This should be your column header in your excel, should be exactly to be the same as in your excel, otherwise it won't read in these covariates.
    name{2}='';
    val=[];
    for v=1:length(name)
        val{v}=M{:,name{v}};
    end
end

%%%%%%%%%%%%%%%%%%%%%%should not edit below %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%make output folder for each contrast
if ~exist(out_dir)
    mkdir(out_dir);
end
cd(out_dir);

%covariates
%pass the covariates to a struct
if cov==1
    covariates.name=name;
    for i=1:cov_num
        values{i}=transpose(val{i});
    end
    covariates.values=values;
else
    covariates={};
end


if test==1 || test==2
    %load the contrast file path for each subject
    scan=[];
    contrast=[];
    for i=1:length(subjects)
        deweight_spm=[root '/' subjects{i} '/' analysis_folder '/SPM.mat'];
        deweight_p=fileparts(deweight_spm);
        load(deweight_spm);
        contrast_names=[];
        scan_files=[];
        for ii=1:length(SPM.xCon)
            contrast_names{ii,1}=SPM.xCon(ii).name;
            scan_files{ii,1}=[deweight_p '/' SPM.xCon(ii).Vcon.fname];
        end
        contrast{i}=contrast_names;
        scan{i}=scan_files;
    end

for ii=1:length(contrast{1})
    out_dirs{ii}=[out_dir '/' contrast{1}{ii}];
    if ~exist(out_dirs{ii})
        mkdir(out_dirs{ii});
    end
end

    allscans=[];
    for i=1:length(scan{1})
        for j=1:length(subjects)
            allscans{i}{j,1}=[scan{j}{i} ',1'];
        end
    end

    if test==1 % one-sample t test
        onesample_t(out_dirs,allscans,covariates);

    elseif test==2 %multiple regression analysis
        multiple_regression(out_dirs,allscans,covariates);

    end

elseif test==3
    for m=1:length(subjects1)
        deweight_spm=[root '/' subjects1{m} '/' analysis_folder '/SPM.mat'];
        deweight_p=fileparts(deweight_spm);
        load(deweight_spm);
        contrast_names1=[];
        scan_files1=[];
        for ii=1:length(SPM.xCon)
            contrast_names1{ii,1}=SPM.xCon(ii).name;
            scan_files1{ii,1}=[deweight_p '/' SPM.xCon(ii).Vcon.fname];
        end
        contrast1{m}=contrast_names1;
        scan1{m}=scan_files1;
    end

for ii=1:length(contrast1{1})
    out_dirs{ii}=[out_dir '/' contrast1{1}{ii}];
    if ~exist(out_dirs{ii})
        mkdir(out_dirs{ii});
    end
end    
        scans1=[];
    for i=1:length(scan1{1})
        for j=1:length(subjects1)
            scans1{i}{j,1}=[scan1{j}{i} ',1'];
        end
    end

    for m=1:length(subjects2)
        deweight_spm=[root '/' subjects2{m} '/' analysis_folder '/SPM.mat'];
        deweight_p=fileparts(deweight_spm);
        load(deweight_spm);
        contrast_names2=[];
        scan_files2=[];
        for ii=1:length(SPM.xCon)
            contrast_names2{ii,1}=SPM.xCon(ii).name;
            scan_files2{ii,1}=[deweight_p '/' SPM.xCon(ii).Vcon.fname];
        end
        contrast2{m}=contrast_names2;
        scan2{m}=scan_files2;
    end
  
        scans2=[];
    for i=1:length(scan2{1})
        for j=1:length(subjects2)
            scans2{i}{j,1}=[scan2{j}{i} ',1'];
        end
    end
   twosample_t(out_dirs,scans1,scans2,covariates); 
end





