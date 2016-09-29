clc, clear;

cwd=pwd;
gname = spm_select(1,'mat','Select lookup file','',fullfile(cwd,'Lookups'));
load(gname)
gname = spm_str_manip(gname,'rt');

load(fullfile(pwd,'matfiles','CheckVBM.mat'));

Q = find(hasRS);

ANOVAdir = spm_select(1,'dir','Select DIR of RS-Nifits for single Seed','',fullfile(cwd,'Projects','RS','GLM'));


for i=1:numel(Q)
    matlabbatch{1}.spm.tools.vbm8.tools.check_cov.data{i} = ...
        fullfile(ANOVAdir,[sub{Q(i)} '.nii']);
end
matlabbatch{1}.spm.tools.vbm8.tools.check_cov.nuisance(1).c = Cov(Q,2);
matlabbatch{1}.spm.tools.vbm8.tools.check_cov.nuisance(2).c = Cov(Q,1);

matlabbatch{1}.spm.tools.vbm8.tools.check_cov.gap = 1;

spm_jobman('run',matlabbatch)



