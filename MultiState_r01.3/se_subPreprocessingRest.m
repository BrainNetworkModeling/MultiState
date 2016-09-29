function se_subPreprocessingRest(subject,xsub,FWHM,conversion,TR)

ms_dir = pwd;
dcm2nii = 0;    % set to 1 to skip SPM dcm conversion and use dcm2nii only


B0dir  = fullfile(subject,'DICOM','B0');        %   Fiedmap DICOM directory
dcmdir = fullfile(subject,'DICOM','EPI');       %   Resting-state EPI DICOM folder 
FMdir  = fullfile(subject,'Fieldmap');          %   Fiedmap directory
origdir  = fullfile(subject,'RS','Orig');       %   Unprocessed resting-state EPI folder
normdir  = fullfile(subject,'RS','Normalised'); %   Normalized resting-state EPI folder
ICAdir  = fullfile(subject,'RS','ICA');         %   Smoothed & normalized 4Dzipped EPIs for MELODIC

% --------------- Convert source files ---------------------

if conversion == 0 
     % ---------- split 4D NIFTI ---------------
    nifti_file = fullfile(dcmdir,[xsub '.nii.gz']);
    mkdir(origdir);
    gunzip(nifti_file);
    nifti_file = fullfile(dcmdir,[xsub '.nii']);
    spm_file_split(nifti_file,origdir);
    delete(nifti_file);
    % ---------- delete dummys --------------------
    fils = dir(fullfile(origdir,'*.nii')); fils = fils(~[fils.isdir]);
    for ex = 1:4
        delete(fullfile(origdir,fils(ex).name))
    end
   
    if exist(B0dir) && ~exist(FMdir)
    % ----- Fieldmap -------
        nifti_gz = dir(fullfile(B0dir,'*.nii.gz'));
        mkdir(FMdir);
        for i=1:numel(nifti_gz)
            gunzip(fullfile(B0dir,nifti_gz(i).name));
        end
        nifti_file = dir(fullfile(B0dir,'*.nii'));
        for i=1:numel(nifti_file)
            spm_file_split(fullfile(B0dir,nifti_file(i).name),FMdir);
            delete(fullfile(B0dir,nifti_file(i).name));
        end
    end


elseif conversion == 1
    % ---------- DICOM conversion --------------------
    mkdir(origdir);
 
 % -----    SPM conversion    -----
    try; if dcm2nii == 0
        load(fullfile(ms_dir,'matfiles','se_Convert.mat')); matlabbatch{1}.spm.util.dicom.data = {};
        fils = dir(fullfile(dcmdir,'*')); fils = fils(~[fils.isdir]);
        if numel(fils)<176
            fils = dir(fullfile(dcmdir,'*.*')); fils = fils(~[fils.isdir]);
        end
        for file = 1:size(fils,1)
            matlabbatch{1}.spm.util.dicom.data{file,1} = fullfile(dcmdir,fils(file).name);
        end
        matlabbatch{1}.spm.util.dicom.outdir{1} = origdir;
        spm_jobman('run',matlabbatch)
    end; end

 % --------- dcm2nii conversion by Chris Rorden :: 6 June 2013 -----    
    if numel(dir(fullfile(origdir,'*.nii'))) <= 150    
        if numel(strfind(computer,'MAC'))>0
            system([fullfile(ms_dir,'matfiles','dcm2niiMAC') ' -b ' fullfile(ms_dir,'matfiles','dcm2nii.ini') ...
                  ' -o ' origdir ' ' fullfile(dcmdir,'*.*')]);
        else
            system([fullfile(ms_dir,'matfiles','dcm2nii') ' -b ' fullfile(ms_dir,'matfiles','dcm2nii.ini') ...
                  ' -o ' origdir ' ' fullfile(dcmdir,'*.*')]);
        end
    end
% ---------- delete dummys --------------------
    fils = dir(fullfile(origdir,'*.nii')); fils = fils(~[fils.isdir]);
    for ex = 1:4
        delete(fullfile(origdir,fils(ex).name))
    end


    if exist(B0dir) && ~exist(FMdir)
% ----- Fieldmap -----
        mkdir(FMdir);
        if numel(strfind(computer,'MAC'))>0
            system([fullfile(ms_dir,'matfiles','dcm2niiMAC') ' -b ' fullfile(ms_dir,'matfiles','dcm2nii.ini') ...
                  ' -o ' FMdir ' ' fullfile(B0dir,'*.*')]);
        else
            system([fullfile(ms_dir,'matfiles','dcm2nii') ' -b ' fullfile(ms_dir,'matfiles','dcm2nii.ini') ...
                  ' -o ' FMdir ' ' fullfile(B0dir,'*.*')]);
        end
    end
end
    

%%% -------- Proceed with pre-processing -------------------


if numel(dir(fullfile(subject,'RS',[xsub '_meanEPI.nii'])))==0

% ---------- realignment (& unwarp) --------------------

    if exist(FMdir)
        if isempty(dir(fullfile(FMdir,'vdm*')))
    % ----- create field map (vdm file) -----
            load(fullfile(ms_dir,'matfiles','UKAfieldmap.mat'))
            fils = dir(fullfile(FMdir,'*.nii')); fils = fils(~[fils.isdir]);
            matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.phase{1} = fullfile(FMdir,fils(3).name);
            matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.magnitude{1} = fullfile(FMdir,fils(1).name);
            matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.mflags.template{1} = ...
                deblank(fullfile(spm('dir'),'templates','T1.nii'));
            fils = dir(origdir); fils = fils(~[fils.isdir]);
            matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.session.epi{1} = fullfile(origdir,fils(1).name);
            matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.anat{1} = fullfile(subject,'3D',[xsub '.nii']);
            spm_jobman('run',matlabbatch)
            movefile(fullfile(origdir,'u*.nii'),fullfile(subject,'RS'));  % delete unwarped Niftis
            movefile(fullfile(origdir,'w*.nii'),fullfile(subject,'RS'));  % delete unwarped Niftis
        else
            fprintf('\n Fieldmap already calculated; skipping calculation...\n')
        end        
            
    % ----- realign & unwarp (using the vdm file) -----
        try
            load(fullfile(ms_dir,'matfiles','realign_unwarp.mat'))
            matlabbatch{1}.spm.spatial.realignunwarp.data.scans = {};
            fils = dir(fullfile(origdir,'*.*')); fils = fils(~[fils.isdir]);
            for file = 1:size(fils,1)
                matlabbatch{1}.spm.spatial.realignunwarp.data.scans{file,1} = fullfile(origdir,fils(file).name);
            end
            Fieldmap = dir(fullfile(FMdir,'vdm5*.nii'));
            matlabbatch{1}.spm.spatial.realignunwarp.data.pmscan{1} = fullfile(FMdir,Fieldmap.name);
            spm_jobman('run',matlabbatch)

            mEPI = dir(fullfile(origdir,'mean*.nii'));
            movefile(fullfile(origdir,mEPI.name),...
                fullfile(subject,'RS',[xsub '_meanEPI.nii']))
            copyfile(fullfile(origdir,'*.txt'),...
                fullfile(subject,'RS'))
            delete(fullfile(origdir,[xsub '_*.nii']));  % delete unwarped Niftis
        catch
            fprintf('Problem with Fieldmap of %s: %s \n', xsub, err.message);
        end
    end

    if numel(dir(fullfile(subject,'RS',[xsub '_meanEPI.nii'])))==0
    % ---------- realignment --------------------
        load(fullfile(ms_dir,'matfiles','se_realign.mat'))
        matlabbatch{1}.spm.spatial.realign.estwrite.data = {};
        fils = dir(fullfile(origdir,'*.*')); fils = fils(~[fils.isdir]);
        for file = 1:size(fils,1)
            matlabbatch{1}.spm.spatial.realign.estwrite.data{1}{file,1} = fullfile(origdir,fils(file).name);
        end
        spm_jobman('run',matlabbatch)
        
        
        mEPI = dir(fullfile(origdir,'mean*.nii'));
        movefile(fullfile(origdir,mEPI.name),...
            fullfile(subject,'RS',[xsub '_meanEPI.nii']))
        copyfile(fullfile(origdir,'*.txt'),...
            fullfile(subject,'RS'))
    end
end


if numel(dir(fullfile(subject,'RS',[xsub '_meanEPI_seg_sn.mat'])))==0
    
    % ---------- Affine registration into MNI space --------------------
    load(fullfile(ms_dir,'matfiles','se_coregisterToMNI.mat')); cnt = 1;
    matlabbatch{1}.spm.spatial.coreg.estimate.source{1} = fullfile(subject,'RS',[xsub '_meanEPI.nii']);
    matlabbatch{1}.spm.spatial.coreg.estimate.ref{1}    = fullfile(spm('dir'),'tpm','grey.nii');
    fils = dir(fullfile(origdir,'*.nii')); fils = fils(~[fils.isdir]); cnt = 1;
    for file = 1:size(fils,1)
        matlabbatch{1}.spm.spatial.coreg.estimate.other{cnt,1} = fullfile(origdir,fils(file).name);  cnt = cnt+1;
    end
    spm_jobman('run',matlabbatch)

    
    % ---------- Normalise using unified segmentation -----
    load(fullfile(ms_dir,'matfiles','se_segment.mat'));
    matlabbatch{1}.spm.spatial.preproc.data{1} = fullfile(subject,'RS',[xsub '_meanEPI.nii']);
    matlabbatch{1}.spm.spatial.preproc.opts.tpm{1} = deblank(fullfile(spm('dir'),'tpm','grey.nii'));
    matlabbatch{1}.spm.spatial.preproc.opts.tpm{2} = deblank(fullfile(spm('dir'),'tpm','white.nii'));
    matlabbatch{1}.spm.spatial.preproc.opts.tpm{3} = deblank(fullfile(spm('dir'),'tpm','csf.nii'));
    spm_jobman('run',matlabbatch)
end


if numel(dir(normdir))==0 && numel(dir(fullfile(subject,'RS','Orig.zip')))==0
    
 % ---------- Apply deformations --------------------
    load(fullfile(ms_dir,'matfiles','se_normalise.mat'));
    mkdir(normdir)
    matlabbatch{1}.spm.util.defs.comp{1}.sn2def.matname{1} = fullfile(subject,'RS',[xsub '_meanEPI_seg_sn.mat']);
    matlabbatch{1}.spm.util.defs.comp{2}.sn2def.matname{1} = fullfile(ms_dir,'matfiles','MNI152_T1_2mm_seg_inv_sn.mat');
    matlabbatch{1}.spm.util.defs.comp{2}.sn2def.vox = [NaN NaN NaN];    
    
    fils = dir(fullfile(origdir,'*.nii*')); fils = fils(~[fils.isdir]);
    matlabbatch{1}.spm.util.defs.fnames = {};
    for file = 1:size(fils,1)
        matlabbatch{1}.spm.util.defs.fnames{file,1} = fullfile(origdir,fils(file).name);
    end
    matlabbatch{1}.spm.util.defs.savedir.saveusr{1}        = normdir;
    matlabbatch{1}.spm.util.defs.interp = 3;
    spm_jobman('run',matlabbatch)
    
    matlabbatch{1}.spm.util.defs.fnames = {};
    matlabbatch{1}.spm.util.defs.fnames{1} =  fullfile(subject,'RS',[xsub '_meanEPI.nii']);
    matlabbatch{1}.spm.util.defs.savedir.saveusr{1}        = fullfile(subject,'RS');
    matlabbatch{1}.spm.util.defs.interp = 1;
    spm_jobman('run',matlabbatch)
    
end


if numel(dir(normdir))~=0

 % ---------- Smooth --------------------
    load(fullfile(ms_dir,'matfiles','se_smooth.mat'));
    mkdir(fullfile(subject,'RS',['Smooth_' int2str(FWHM) 'mm']));
    matlabbatch{1}.spm.spatial.smooth.fwhm = [FWHM FWHM FWHM];
    fils = dir(fullfile(normdir,'w*.nii*')); fils = fils(~[fils.isdir]);
    for file = 1:size(fils,1)
        matlabbatch{1}.spm.spatial.smooth.data{file,1} = fullfile(normdir,fils(file).name);
    end
    spm_jobman('run',matlabbatch)
    movefile(fullfile(normdir,'sw*.nii*'),fullfile(subject,'RS',['Smooth_' int2str(FWHM) 'mm']))
 
    
    
 % ---------- cleanup --------------------
    zip(fullfile(subject,'RS',['Orig.zip']),origdir)
    [~, ~ ] = rmdir(normdir,'s');

end
    
    
if numel(dir(fullfile(subject,'RS',['c1w' xsub '_meanEPI.nii'])))==0
   
 % ---------- Get final segments from T1 & EPI --------------------
    load(fullfile(ms_dir,'matfiles','se_NewSegment.mat'))
    for i=1:6;
        matlabbatch{1}.spm.tools.preproc8.tissue(i).tpm{1} = fullfile(spm('dir'),'toolbox','Seg',['TPM.nii,' int2str(i)]);
    end
    matlabbatch{1}.spm.tools.preproc8.channel.vols    = {};
    matlabbatch{1}.spm.tools.preproc8.channel.vols{1} = fullfile(subject,'RS',['w' xsub '_meanEPI.nii']);
    if numel(dir(fullfile(subject,'3D',['wmr' xsub '.nii'])))==1
        matlabbatch{1}.spm.tools.preproc8.channel.vols{2} = fullfile(subject,'3D',['wmr' xsub '.nii']);
    end
    spm_jobman('run',matlabbatch)
    delete(fullfile(subject,'RS','w*_seg8.mat'))

end


if numel(dir(fullfile(subject,'RS',['MNI152GM_' xsub '.nii'])))==0 || ...
        numel(dir(fullfile(subject,'RS',['Counfounds_' xsub '.mat'])))==0
    clear pXYZ tX tY tZ gx gx2 rp rpp A dat msk gx2 reg y

    for i=1:3
        msk{i} = spm_read_vols(spm_vol(fullfile(subject,'RS',['c' int2str(i) 'w' xsub '_meanEPI.nii'])));
    end
    msk{4} = (msk{1}+msk{2}+msk{3});
    
    fil = dir(fullfile(subject,'RS',['Smooth_' int2str(FWHM) 'mm'],'*.nii'));
    
    try
    txt  = dir(fullfile(subject,'RS','*.txt'));
        rp   = load(fullfile(subject,'RS',txt(1).name));
    catch
        txt  = dir(fullfile(origdir,'*.txt'));
        rp   = load(fullfile(origdir,txt(1).name));
    end
    rp   = rp-repmat(nanmean(rp),size(rp,1),1);
    rpp  = [zeros(1,6); diff(rp)];

    
    for i=1:numel(fil); Vi(i)  = spm_vol(fullfile(subject,'RS',['Smooth_' int2str(FWHM) 'mm'],fil(i).name)); end

    gx     = nan(3,numel(fil));
    ind    = find((msk{4})>.5);
    y      = zeros(numel(ind),numel(Vi));
    
    for i=1:numel(fil)
        dat    = spm_read_vols(Vi(i));
        y(:,i) = dat(ind);
        for ii=1:4
            A = dat.*msk{ii}; A = A(msk{ii}>.1); gx(ii,i) = sum(A)/sum(msk{ii}(msk{ii}>.1));
        end
    end
    
    y = y - repmat(mean(y),size(y,1),1);
    COEFF = princomp(y,'econ');

    gx2(1,:) = gx(1,:)-mean(gx(1,:));
    gx2(2,:) = gx(2,:)-mean(gx(2,:));
    gx2(3,:) = gx(3,:)-mean(gx(3,:));
    gx2(4,:) = gx(4,:)-mean(gx(4,:));
    
    reg = [gx2' gx2'.^2 rp rp.^2 rpp rpp.^2 COEFF(:,1:5)];
    reg = reg-repmat(nanmean(reg),size(reg,1),1);
    reg = reg./repmat(nanstd(reg),size(reg,1),1);
    
    if conversion == 1
        try
            DICOMS = [dir(fullfile(dcmdir,'*.*'));dir(fullfile(dcmdir,'*'))]; DICOMS = DICOMS(~[DICOMS.isdir]);
            TR = dicominfo(fullfile(dcmdir,DICOMS(1).name));
            TR = TR.RepetitionTime/1000;
        end
    end
    
    save(fullfile(subject,'RS',['Counfounds_' xsub '.mat']),'Vi','reg','rp','rpp','gx','gx2','TR')
    
    
    load(fullfile(ms_dir, 'MaskenEtc/MNI152GM.mat'));
    XYZ = Vi(1).mat \ VM.mat * [maskXYZ; ones(1,size(maskXYZ,2))];
    dat = spm_get_data(Vi,XYZ); clear V
    
    V = struct('fname',	fullfile(subject,'RS',['MNI152GM_' xsub '.nii']),...
        'dim',		[numel(Vi) size(XYZ,2) 1],...
        'dt',       [8 0],...
        'mat',		eye(4),...
        'descrip',	[xsub ' 5 mm']);
    
    V       = spm_write_vol(V,dat);
    save(fullfile(subject,'RS',['Counfounds_' xsub '.mat']),'Vi','V','reg','rp','rpp','gx','gx2','TR')
end


try; rmdir(fullfile(subject,'RS',['Smooth_' int2str(FWHM) 'mm']),'s'); end



if numel(dir(fullfile(ICAdir,[xsub '.nii.gz'])))==0
    
    n = spm_vol(fullfile(subject,'RS',['MNI152GM_' xsub '.nii'])); n = n.dim(1);
    [~, ~ ] = rmdir(ICAdir,'s');
    fils = dir(fullfile(origdir,'*.nii*'));
    
 % ---------- Apply deformations --------------------
    load(fullfile(ms_dir,'matfiles','se_normalise.mat'));
    mkdir(ICAdir)
    matlabbatch{1}.spm.util.defs.comp{1}.sn2def.matname{1} = fullfile(subject,'RS',[xsub '_meanEPI_seg_sn.mat']);
    matlabbatch{1}.spm.util.defs.comp{2}.sn2def.matname{1} = fullfile(ms_dir,'matfiles','MNI152_T1_2mm_seg_inv_sn.mat');
    matlabbatch{1}.spm.util.defs.comp{2}.sn2def.vox = [3 3 3];    
    
    fils = dir(fullfile(origdir,'*.nii*')); fils = fils(~[fils.isdir]);
    try
        fils = fils(end+1-n:end);
    catch
        fils = fils(end-n:end);
    end
    matlabbatch{1}.spm.util.defs.fnames = {};
    for file = 1:size(fils,1)
        matlabbatch{1}.spm.util.defs.fnames{file,1} = fullfile(origdir,fils(file).name);
    end
    matlabbatch{1}.spm.util.defs.savedir.saveusr{1}        = ICAdir;
    matlabbatch{1}.spm.util.defs.interp = 3;
    spm_jobman('run',matlabbatch)

    
    
 % ---------- Smooth --------------------
    load(fullfile(ms_dir,'matfiles','se_smooth.mat'));
    matlabbatch{1}.spm.spatial.smooth.fwhm = [FWHM FWHM FWHM];
    fils = dir(fullfile(ICAdir,'w*.nii*')); fils = fils(~[fils.isdir]);
    for file = 1:size(fils,1)
        matlabbatch{1}.spm.spatial.smooth.data{file,1} = fullfile(ICAdir,fils(file).name);
    end
    spm_jobman('run',matlabbatch)
    for file = 1:size(fils,1)
        delete(fullfile(ICAdir,fils(file).name));
    end
    
    if numel(strfind(computer,'MAC'))>0
        setenv('FSLDIR','/usr/local/fsl');
        setenv('FSLOUTPUTTYPE', 'NIFTI_GZ');
    end

    system([getenv('FSLDIR') '/bin/fslmerge -t ' fullfile(ICAdir,[xsub '.nii']) ' ' fullfile(ICAdir,'sw*.nii')]);

 % ---------- cleanup --------------------
    rmdir(origdir,'s')
    delete(fullfile(ICAdir,'sw*.nii'));
end
    
    
