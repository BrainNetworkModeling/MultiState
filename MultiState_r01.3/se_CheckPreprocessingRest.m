clear, clc; warning off all; spm('Defaults','fmri')


force = 0;  % Try to preprocess all subjects, regardless of previous preprocessing

%--------------- readout of path-to-data ----------------------------
fid = fopen('_path2data','r'); datadir = fgetl(fid); fclose(fid);


Di = spm_select(Inf,'dir','Select subject directories','',fullfile(datadir,'DATA'));

for i=1:size(Di,1)
    subjects{i} = strrep(Di(i,:),' ','');
    nx = strrep(Di(i,:),' ','');
    if nx(end) == filesep; nx(end) = []; end
    us = strfind(nx,filesep);
    nx = nx(us(end)+1:end);
    xsub{i} = nx;
end


TR = 0.72; % input('Enter TR of cohort in sec (Ex.: 2.2): ');


% First determine which subjects need processing, to ensure an even
% distribution across workers

idx=[];
i=1;

if force
    idx = 1 : numel(subjects);
else
    for sub = 1:numel(subjects)
        if numel(dir(fullfile(subjects{sub},'RS','ICA',[xsub{sub} '.nii.gz'])))==0
            idx(i) = sub;
            i = i + 1;
        end
    end
end

fprintf('Found %i unprocessed of %i subjects..\n', numel(idx), numel(subjects));

if numel(idx) == 0
    fprintf('No subjects to process! Aborting.');
    return; 
end
subjects = subjects(idx);
xsub = xsub(idx);


parfor sub = 1:numel(subjects)
% What conversion is necessary? 0 = unzip & split 4D NIFTI, 1 = DICOM or PARREC conversion here
    conversion = -1;
    fils = [];
    
    subdcmdir = fullfile(subjects{sub},'DICOM','EPI');
    niigzsub = fullfile(subdcmdir,[xsub{sub} '.nii.gz']);
    niftis = fullfile(subjects{sub},'RS','Orig','*.nii');
    
    niigzstart = dir(niigzsub); 
    niftistart = dir(niftis);

    
    niigznames = dir(fullfile(subdcmdir,'*.nii.gz'));
    if numel(niigzstart)==0 && numel(niftistart) <= 150 && numel(niigznames)>=1
    % rename zipped 4D nifti file
        movefile(fullfile(subdcmdir,niigznames(1).name),niigzsub);
        conversion = 0;
    end

    if conversion == -1 && numel(dir(fullfile(subdcmdir,'*.nii'))) >= 1
        if numel(dir(fullfile(subdcmdir,'*.nii'))) >= 150
    % merge & zip nifti file
            system([getenv('FSLDIR') '/bin/fslmerge -t ' niigzsub ' ' fullfile(subdcmdir,'*.nii')]);
        elseif numel(dir(fullfile(subdcmdir,'*.nii'))) == 1
    % zip nifti file
            niftiname= dir(fullfile(subdcmdir,'*.nii'));
            gzip(fullfile(subdcmdir,niftiname.name));
            try movefile(fullfile(subdcmdir,[niftiname.name '.gz']),niigzsub); end
        end
        if numel(dir(niigzsub))==1
            delete(fullfile(subdcmdir,'*.nii'));
            conversion = 0;
        end
    end

    if conversion == -1 && numel(dir(niigzsub))==1 && numel(dir(niftis)) <= 150
        conversion = 0;
    elseif conversion == -1 && numel(dir(niftis)) <= 150 &&...
            numel(dir(fullfile(subdcmdir,[xsub{sub} '_nodummies.nii.gz'])))==1
        conversion = 2;
    elseif conversion == -1 && numel(dir(niftis)) > 150
        conversion = 3;       
    end

    if conversion == -1;
    % Try DICOM (PARREC)
        fils = dir(fullfile(subdcmdir,'*')); fils = fils(~[fils.isdir]);
        if numel(fils)>150
            conversion = 1;
        else
            fils = dir(fullfile(subdcmdir,'*.*')); fils = fils(~[fils.isdir]);
        end
        if numel(fils)>150
            conversion = 1;
        elseif numel(dir(fullfile(subdcmdir,'*.PAR')))==1 ||...
                numel(dir(fullfile(subdcmdir,'*.par')))==1
            conversion = 1;
        else
            fprintf('No viable images found for subject %s. Skipping...\n', xsub{sub});
        end
    end
    
    % ----- Start preprocessing -----  
    try
        if conversion >= 0
            se_subPreprocessingRest(subjects{sub},xsub{sub},5,conversion,TR);
        else
            fprintf('No viable images found for subject %s. Skipping...\n', xsub{sub});
        end
    catch err
        fprintf('Error processing %s: %s\n', xsub{sub}, err.message);
    end
end
    
    
