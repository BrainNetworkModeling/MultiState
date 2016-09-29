function se_CheckPreprocessingRest_nogui ( list_file, TR, force )

% Accepts a list of subjects for preprocessing of resting state fMRI.
% 
% list_file      - path to text file containing a list of subjects to process
% TR             - repetition time, in seconds
% force          - whether to force reprocessing of existing subjects
%                  [default = false]

%fprintf('Current dir: %s\n',pwd);

if nargin < 3
	force = 0;
end

if ~exist('TR','var')
	error('Please specify a TR!\n');
end

fprintf('Loading subjects from %s\n', list_file);

warning off all; spm('Defaults','fmri')

% Read subjects from text file
Di = textread(list_file, '%s', 'delimiter', '\n');
Di = char(Di);

for i=1:size(Di,1)
    subjects{i} = strrep(Di(i,:),' ','');
    nx = strrep(Di(i,:),' ','');
    if nx(end) == filesep; nx(end) = []; end
    us = strfind(nx,filesep);
    nx = nx(us(end)+1:end);
    xsub{i} = nx;
end

% First determine which subjects need processing, to ensure an even
% distribution across workers

clear idx;
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
    
    subdcmdir = fullfile(subjects{sub},'DICOM','EPI');      %   Subjects resting-state EPI DICOM folder 
    niigzsub = fullfile(subdcmdir,[xsub{sub} '.nii.gz']);   %   Subjects 4Dzipped EPI file
    niftis = fullfile(subjects{sub},'RS','Orig','*.nii');   %   Subjects 3D EPI files
     
    niigzstart = dir(niigzsub); 
    niftistart = dir(niftis);

    
    niigznames = dir(fullfile(subdcmdir,'*.nii.gz'));
    if numel(niigzstart)==0 && numel(niftistart) <= 170 && numel(niigznames)>=1
    % rename zipped 4D nifti file
        movefile(fullfile(subdcmdir,niigznames(1).name),niigzsub);
        conversion = 0;
    end

    if conversion == -1 && numel(dir(fullfile(subdcmdir,'*.nii'))) >= 1
        if numel(dir(fullfile(subdcmdir,'*.nii'))) >= 170
        % merge & zip nifti file
            system([getenv('FSLDIR') '/bin/fslmerge -t ' niigzsub ' ' fullfile(subdcmdir,'*.nii')]);
        elseif numel(dir(fullfile(subdcmdir,'*.nii'))) == 1
            niftiname= dir(fullfile(subdcmdir,'*.nii'));
            gzip(fullfile(subdcmdir,niftiname.name));
            try movefile(fullfile(subdcmdir,[niftiname.name '.gz']),niigzsub); end
        end
        if numel(dir(niigzsub))==1
            delete(fullfile(subdcmdir,'*.nii'));
            conversion = 0;
        end
    end

    if conversion == -1 && numel(dir(niigzsub))==1 && numel(dir(niftis)) <= 170
        conversion = 0;
    elseif conversion == -1 && numel(dir(niftis)) <= 170 &&...
            numel(dir(fullfile(subdcmdir,[xsub{sub} '_nodummies.nii.gz'])))==1
        conversion = 2;
    elseif conversion == -1 && numel(dir(niftis)) > 170
        conversion = 3;       
    end

    if conversion == -1;
    % Try DICOM (PARREC)
        fils = dir(fullfile(subdcmdir,'*.*')); fils = fils(~[fils.isdir]);
        if numel(fils)>170
            conversion = 1;
        else
            fils = dir(fullfile(subdcmdir,'*.*')); fils = fils(~[fils.isdir]);
        end
        if numel(fils)>170
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
