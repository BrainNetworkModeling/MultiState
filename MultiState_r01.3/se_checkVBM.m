clc, clear;

cwd=pwd;
gname = spm_select(1,'mat','Select lookup file','',fullfile(cwd,'Lookups'));
load(gname)
gname = spm_str_manip(gname,'rt');


load(fullfile(pwd,'matfiles','CheckVBM.mat'));

c = spm_input('Which images',1,'m','sm0|m0|sm|m',[1 2 3 4],1);
C = {'sm0','m0','sm','m'};

pref = C{c};

Q = find(hasVBM);

for i=1:numel(Q)
    matlabbatch{1}.spm.tools.vbm8.tools.check_cov.data{i} = ...
        fullfile(pwd,'DATA',SubDir{Q(i)},sub{Q(i)},'3D',[pref 'wrp1'  sub{Q(i)} '.nii']);
end
matlabbatch{1}.spm.tools.vbm8.tools.check_cov.nuisance(1).c = Cov(Q,2);
matlabbatch{1}.spm.tools.vbm8.tools.check_cov.nuisance(2).c = Cov(Q,1);

matlabbatch{1}.spm.tools.vbm8.tools.check_cov.gap = 1;
if c>2
    matlabbatch{1}.spm.tools.vbm8.tools.check_cov.scale = 1;
end

spm_jobman('run',matlabbatch)



