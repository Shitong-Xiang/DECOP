clear;clc

path = '/share/inspurStorage/home1/ISTBI_data/ABCD_MID/First_level_new';
addpath(genpath('/share/inspurStorage/home2/jiaozy/ABCD_Task_pipline/spm12'));

files = dir(fullfile(path,'sub-*'));
error = zeros(1,length(files));
parfor i = 1:length(files)
    disp(i);
    con = dir(fullfile(path,files(i).name,'con*'));
    if isempty(con)
        error(i) = i;
    end
end
error = error(error ~= 0);
files(error) = [];

SubID = {files.name};
condition = {'con_0001.nii','con_0005.nii','con_0009.nii','con_0017.nii','con_0013.nii'};

target = '/share/inspurStorage/home1/ISTBI_data/ABCD_MID/FD_QC/Decode';
mkdir(target)
cd(target)

subs = cell(1,length(SubID));
for k= 1:length(SubID)
    disp(['First',num2str(k)])
    name = SubID{k};
    output = fullfile(target,name);
    matlabbatch{1}.spm.stats.factorial_design.dir = cellstr(output);
    
    for j = 1:length(condition)
        matlabbatch{1}.spm.stats.factorial_design.des.mreg.scans(j,1)= {fullfile(path,name,condition{j})};
    end
    
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(1).c = [-2 -1 0 1 2]';
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(1).cname = 'Value';
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(1).iCC = 1;
    
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(2).c = [2 1 0 1 2]';
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(2).cname = 'Salience';
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(2).iCC = 1;
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.incint = 1;
    
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(3).c = [1 -2 2 -2 1]';
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(3).cname = 'W';
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(3).iCC = 1;
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.incint = 1;
    
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(4).c = [-1 2 0 -2 1]';
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(4).cname = 'N';
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(4).iCC = 1;
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.incint = 1;
    
    matlabbatch{1}.spm.stats.factorial_design.masking.em = {'/share/inspurStorage/home1/ISTBI_data/ABCD_MID/MNI152_T1_3mm_brain_mask.nii'};
    matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
    
    
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
        
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'mean';
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.convec = [1 0 0 0 0];
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'Value';
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.convec = [0 1 0 0 0];
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    
    matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = 'Salience';
    matlabbatch{3}.spm.stats.con.consess{3}.tcon.convec = [0 0 1 0 0];
    matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
    
    matlabbatch{3}.spm.stats.con.consess{4}.tcon.name = 'W';
    matlabbatch{3}.spm.stats.con.consess{4}.tcon.convec = [0 0 0 1 0];
    matlabbatch{3}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
    
    matlabbatch{3}.spm.stats.con.consess{5}.tcon.name = 'N';
    matlabbatch{3}.spm.stats.con.consess{5}.tcon.convec = [0 0 0 0 1];
    matlabbatch{3}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
    
    subs{k} = matlabbatch;
end

decode_error = zeros(1,length(files));
parfor k = 1:length(SubID)
    try
        spm_jobman('run',subs{k});
    catch
        decode_error(k) = k;
    end
end
decode_error = decode_error(decode_error ~= 0);



