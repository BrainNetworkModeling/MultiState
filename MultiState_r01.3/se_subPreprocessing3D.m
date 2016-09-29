function se_subPreprocessing3D(subject,xsub,structurals,conversion,force)
ms_dir = pwd;
dcm2nii = 0;

subdcmdir = fullfile(subject,'DICOM',structurals);
sub3Ddir  = fullfile(subject,structurals);

if conversion == 1
    % ---------- DICOM conversion --------------------
    mkdir(sub3Ddir);
    
    % -----    SPM conversion    -----
    try; if dcm2nii == 0
        load(fullfile(ms_dir,'matfiles','se_Convert.mat')); matlabbatch{1}.spm.util.dicom.data = {};    
        fils = dir(fullfile(subdcmdir,'*')); fils = fils(~[fils.isdir]);
        if numel(fils)<100
            fils = dir(fullfile(subdcmdir,'*.*')); fils = fils(~[fils.isdir]);
        end

        if numel(fils)>90
            for file = 1:size(fils,1)
                matlabbatch{1}.spm.util.dicom.data{file,1} = fullfile(subdcmdir,fils(file).name);
            end
            matlabbatch{1}.spm.util.dicom.outdir{1} = sub3Ddir;
            spm_jobman('run',matlabbatch)

            fx = dir(fullfile(sub3Ddir,'*.nii'));
            movefile(fullfile(sub3Ddir,fx(find([fx.bytes]==max([fx.bytes]))).name),fullfile(sub3Ddir,[xsub '.nii']))
            try delete(fx(find([fx.bytes]~=max([fx.bytes]))).name); end
        end
    end; end

    % --------- dcm2nii conversion by Chris Rorden :: 6 June 2013 -----    
    if numel(dir(fullfile(sub3Ddir,'*.nii'))) <= 1   
        if numel(strfind(computer,'MAC'))>0
            system([fullfile(ms_dir,'matfiles','dcm2niiMAC') ' -b ' fullfile(ms_dir,'matfiles','dcm2nii.ini') ...
                  ' -o ' sub3Ddir ' ' fullfile(subdcmdir,'*.*')]);
        else
            system([fullfile(ms_dir,'matfiles','dcm2nii') ' -b ' fullfile(ms_dir,'matfiles','dcm2nii.ini') ...
                  ' -o ' sub3Ddir ' ' fullfile(subdcmdir,'*.*')]);
            fx = dir(fullfile(sub3Ddir,'*.nii'));
            movefile(fullfile(sub3Ddir,fx(find([fx.bytes]==max([fx.bytes]))).name),fullfile(sub3Ddir,[xsub '.nii']))
        end
    end
end

% ------- Onward with pre-processing ----------

if numel(dir(fullfile(sub3Ddir,'rp*.nii')))~=2 || numel(dir(fullfile(sub3Ddir,[xsub '_seg8.mat'])))~=1 || force == 1

    % ---------- Affine registration to IXI-template for robustness --------------------
    load(fullfile(ms_dir,'matfiles','se_coregisterToMNI.mat'));
    matlabbatch{1}.spm.spatial.coreg.estimate.source{1} = fullfile(sub3Ddir,[xsub '.nii']);
    matlabbatch{1}.spm.spatial.coreg.estimate.ref{1}    = ...
        fullfile(spm('dir'),'toolbox','vbm8','avgT1_Dartel_IXI550_MNI152.nii');
    matlabbatch{1}.spm.spatial.coreg.estimate.other = {''};
    spm_jobman('run',matlabbatch)
    
    % ---------- high-dimensional DARTEL normalization --------------------
    load(fullfile(ms_dir,'matfiles','se_VBM8.mat'));
    matlabbatch{1}.spm.tools.vbm8.estwrite.opts.tpm{1} = fullfile(spm('dir'),'toolbox','Seg','TPM.nii');
    matlabbatch{1}.spm.tools.vbm8.estwrite.extopts.dartelwarp.normhigh.darteltpm{1} = ...
        fullfile(spm('dir'),'toolbox','vbm8','Template_1_IXI550_MNI152.nii');
    matlabbatch{1}.spm.tools.vbm8.estwrite.data = {};
    matlabbatch{1}.spm.tools.vbm8.estwrite.data{1,1} = fullfile(sub3Ddir,[xsub '.nii']);
    % matlabbatch{1}.spm.tools.vbm8.estwrite.opts.biasfwhm = 100;     % if intensity is too uniform
    spm_jobman('run',matlabbatch)
end


fx = dir(fullfile(sub3Ddir,'mwrp*.nii'));
if numel(fx)~=2 || diff([fx.bytes]) ~= 0
    load(fullfile(ms_dir,'matfiles','se_VBM8m0.mat'));
    matlabbatch{1}.spm.tools.vbm8.write.output.GM.modulated = 1;
    matlabbatch{1}.spm.tools.vbm8.write.output.WM.modulated = 1;
    matlabbatch{1}.spm.tools.vbm8.write.extopts.dartelwarp.normhigh.darteltpm{1} = ...
        fullfile(spm('dir'),'toolbox','vbm8','Template_1_IXI550_MNI152.nii');
    matlabbatch{1}.spm.tools.vbm8.write.data = {};
    matlabbatch{1}.spm.tools.vbm8.write.data{1,1} = fullfile(sub3Ddir,[xsub '.nii']);
    spm_jobman('run',matlabbatch)
end


fx = dir(fullfile(sub3Ddir,'m0wrp*.nii'));
if numel(fx)~=2 || diff([fx.bytes]) ~= 0 || force == 1
    load(fullfile(ms_dir,'matfiles','se_VBM8m0.mat'));
    matlabbatch{1}.spm.tools.vbm8.write.extopts.dartelwarp.normhigh.darteltpm{1} = ...
        fullfile(spm('dir'),'toolbox','vbm8','Template_1_IXI550_MNI152.nii');
    matlabbatch{1}.spm.tools.vbm8.write.data = {};
    matlabbatch{1}.spm.tools.vbm8.write.data{1,1} = fullfile(sub3Ddir,[xsub '.nii']);
    spm_jobman('run',matlabbatch)
end


fx = dir(fullfile(sub3Ddir,'sm0wrp*.nii'));
if numel(fx)~=2 || diff([fx.bytes]) ~= 0 || force == 1
    load(fullfile(ms_dir,'matfiles','se_smooth.mat'));
    matlabbatch{1}.spm.spatial.smooth.data = {};
    matlabbatch{1}.spm.spatial.smooth.data{1,1} = fullfile(sub3Ddir,['m0wrp1' xsub '.nii']);
    matlabbatch{1}.spm.spatial.smooth.data{2,1} = fullfile(sub3Ddir,['m0wrp2' xsub '.nii']);
    spm_jobman('run',matlabbatch)
end

fx = dir(fullfile(sub3Ddir,'smwrp*.nii'));
if numel(fx)~=2 || diff([fx.bytes]) ~= 0 || force == 1
    load(fullfile(ms_dir,'matfiles','se_smooth.mat'));
    matlabbatch{1}.spm.spatial.smooth.data = {};
    matlabbatch{1}.spm.spatial.smooth.data{1,1} = fullfile(sub3Ddir,['mwrp1' xsub '.nii']);
    matlabbatch{1}.spm.spatial.smooth.data{2,1} = fullfile(sub3Ddir,['mwrp2' xsub '.nii']);
    spm_jobman('run',matlabbatch)
end

delete(fullfile(sub3Ddir,[xsub '_seg8.mat']));

