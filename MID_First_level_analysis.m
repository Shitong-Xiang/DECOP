%% task prepare
clear;clc
path = '/share/inspurStorage/home1/ISTBI_data/ABCD_MID/Prep_finish';
output = '/share/inspurStorage/home1/ISTBI_data/ABCD_MID/First_level_new';
addpath(genpath('/share/inspurStorage/home2/jiaozy/ABCD_Task_pipline/spm12'));

files = dir(fullfile(path,'sub*'));
exist = dir(fullfile(output,'sub*'));
[~,ind,~] = intersect(string({files.name}),string({exist.name}));
files(ind) = [];

fmri_rt = 0.8;
fmri_t = 60;
fmri_t0 = 1;
fmri_Units =  'secs';

wrong_subject = zeros(1,length(files));
sub_job = cell(1,length(files));
for j = 1:length(files)
    try
        disp(['Subject    ',num2str(j,'%012d')]);
       %% condidtions
        onset_file = dir(fullfile(files(j).folder,files(j).name,'func/*onset.mat'));
        clear matlabbatch
        
        matlabbatch{1}.spm.stats.fmri_spec.dir =cellstr(fullfile(output,files(j).name));
        matlabbatch{1}.spm.stats.fmri_spec.timing.units = fmri_Units;
        matlabbatch{1}.spm.stats.fmri_spec.timing.RT = fmri_rt;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = fmri_t;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = fmri_t0;
        
        subpath = fullfile(files(j).folder,files(j).name,'prep/func/');
        motion = dir(fullfile(subpath,'Regress_*.txt'));
        files_nii = dir(fullfile(subpath,'predenoise*nii'));
       
        for i = 1:length(files_nii)
            % nii file
            niifile = spm_select('ExtFPList',files_nii(i).folder,files_nii(i).name);
            matlabbatch{1}.spm.stats.fmri_spec.sess(i).scans = cellstr(niifile);
            % onset
            onset = load(fullfile(onset_file(i).folder,onset_file(i).name));
            conditions = table2struct(onset.onset_all);
            for k = 1:length(conditions)
                matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(k).name = char(conditions(k).type);
                matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(k).tmod = 0;
                matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(k).pmod = struct('name', {}, 'param', {}, 'poly', {});
                matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(k).onset    = conditions(k).onset';
                matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(k).duration = conditions(k).duration';
            end
            % Nuisance Regressors: motion
            matlabbatch{1}.spm.stats.fmri_spec.sess(i).multi = {''};
            matlabbatch{1}.spm.stats.fmri_spec.sess(i).regress = struct('name', {}, 'val', {});
            motionfiles = spm_select('FPList',subpath,[motion(i).name,'.*']); % 6 motion, pre 1s (6) and post 1 s(6), the other 3 is what
            matlabbatch{1}.spm.stats.fmri_spec.sess(i).multi_reg = cellstr(motionfiles);
            matlabbatch{1}.spm.stats.fmri_spec.sess(i).hpf = 128;
            motion_data = importdata(fullfile(subpath,motion(i).name));
            motion_legnth(i) = size(motion_data,2);
        end
        
        % keyboard(Defaults,need't change)
        matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
        matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
        matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
        matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
        matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0;
        matlabbatch{1}.spm.stats.fmri_spec.mask = {'/share/inspurStorage/home1/ISTBI_data/toolbox/Preprocess_pipeline/MNI152/MNI152_T1_3mm_brain.nii'};
        matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
        
        % model estimation
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep;
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tname = 'Select SPM.mat';
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(1).name = 'filter';
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(1).value = 'mat';
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(2).name = 'strtype';
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(2).value = 'e';
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1).sname = 'fMRI model specification: SPM.mat File';
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1).src_output = substruct('.','spmmat');
        matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
         
        % SPM contrasts
        matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep;
        matlabbatch{3}.spm.stats.con.spmmat(1).tname = 'Select SPM.mat';
        matlabbatch{3}.spm.stats.con.spmmat(1).tgt_spec{1}(1).name = 'filter';
        matlabbatch{3}.spm.stats.con.spmmat(1).tgt_spec{1}(1).value = 'mat';
        matlabbatch{3}.spm.stats.con.spmmat(1).tgt_spec{1}(2).name = 'strtype';
        matlabbatch{3}.spm.stats.con.spmmat(1).tgt_spec{1}(2).value = 'e';
        matlabbatch{3}.spm.stats.con.spmmat(1).sname = 'Model estimation: SPM.mat File';
        matlabbatch{3}.spm.stats.con.spmmat(1).src_exbranch = substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1});
        matlabbatch{3}.spm.stats.con.spmmat(1).src_output = substruct('.','spmmat');
       
        % SPM contrast set
        for k = 1:length(conditions)
            if length(motion)> 1
                cond_vector_run1 = zeros(1,20+motion_legnth(1));
                cond_vector_run1(k)= 1;
                cond_vector_run2 = zeros(1,20+motion_legnth(2));
                cond_vector_run2(k)= 1;
                cond_vector_run = [cond_vector_run1,cond_vector_run2];
            else
                cond_vector_run = zeros(1,20+motion_legnth(1));
                cond_vector_run(k)= 1;
            end
            matlabbatch{3}.spm.stats.con.consess{k}.tcon.name = char(conditions(k).type);
            matlabbatch{3}.spm.stats.con.consess{k}.tcon.convec = cond_vector_run;
            matlabbatch{3}.spm.stats.con.consess{k}.tcon.sessrep = 'none';
        end
        sub_job{j}=  matlabbatch;
    catch
        wrong_subject(j) = j;
    end
end
sub_job(wrong_subject~=0)=[];

%% run SPM
wrong_first_subject = zeros(1,length(sub_job));
parfor j = 1:length(sub_job)
    try
        spm_jobman('run',sub_job{j});
    catch
        wrong_first_subject(j) = j;
    end
end
cd(output);
save wrong_sub.mat wrong_subject wrong_first_subject;