function se_CheckPrep3DnRest_nogui( list_file, TR, force )

% Accepts a list of subjects for preprocessing of resting state fMRI.
% 
% list_file      - path to text file containing a list of subjects to process
% force          - whether to force reprocessing of existing subjects
%                  [default = false]

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

%%% ------ VBM PREPROCESSING ------ %%%%%


for i=1:size(Di,1)
    subjects{i} = strrep(Di(i,:),' ','');
    nx = strrep(Di(i,:),' ','');
    if nx(end) == filesep; nx(end) = []; end
    us = strfind(nx,filesep);
    nx = nx(us(end)+1:end);
    xsub{i} = nx;
end

structurals     = {'3D'};



% First determine which subjects need processing, to ensure an even
% distribution across workers

idx=[];
i=1;

if force
    idx = 1 : numel(subjects);
else
    for sub = 1:numel(subjects)
        sub3Ddir  = fullfile(subjects{sub},structurals{1});
        isda = numel(dir(fullfile(sub3Ddir,'sm0wrp*.*'))) ...
                + numel(dir(fullfile(sub3Ddir,'smwrp*.*'))) ...
                + numel(dir(fullfile(sub3Ddir,'m0wrp*.*'))) ...
                + numel(dir(fullfile(sub3Ddir,'mwrp*.*')));
        if isda ~= 8
            idx(i) = sub;
            i = i + 1;
        end

    end
end

fprintf('Found %i unprocessed of %i subjects with 3D data..\n', numel(idx), numel(subjects));

if numel(idx) >= 1

subjects = subjects(idx);
xsub = xsub(idx);


parfor sub = 1:numel(subjects)
    
    % What conversion is necessary? 0 = none (NIFTI exists), 1 = DICOM
    conversion = -1; fils = [];
    
    subdcmdir = fullfile(subjects{sub},'DICOM',structurals{1});     %   T1 DICOM folder
    sub3Ddir  = fullfile(subjects{sub},structurals{1});             %   Subjects VBM folder
    niigzsub = fullfile(subdcmdir,[xsub{sub} '.nii.gz']);           %   Subjects zipped T1 file
    niftisub = fullfile(sub3Ddir,[xsub{sub} '.nii']);               %   Subjects T1 file for VBM
    
    niftistart = dir(niftisub);
    niigzstart = dir(niigzsub); 

    if numel(niftistart)~=1 && numel(niigzstart)~=1
        niftinames= dir(fullfile(subdcmdir,'*.nii'));
        if numel(niftinames) >= 1
    % zip & rename nifti file
            gzip(fullfile(subdcmdir,niftinames(1).name));
            gzniinames= dir(fullfile(subdcmdir,'*.nii.gz'));
            try movefile(fullfile(subdcmdir,gzniinames(1).name),niigzsub); end
            if numel(dir(niigzsub))==1
                delete(fullfile(subdcmdir,niftinames(1).name));        
                conversion = 0;
            end
        end

        niigznames = dir(fullfile(subdcmdir,'*.nii.gz'));
        if conversion == -1 && numel(niigznames)>=1
    % rename zipped nifti file
            movefile(fullfile(subdcmdir,niigznames(1).name),niigzsub);
            conversion = 0;
        end
    end
    
    niigzstart = dir(niigzsub);
    if numel(niigzstart)==1 && numel(niftistart)~=1
    % copy & extract nifti file to 3D dir
        mkdir(fullfile(sub3Ddir));
        copyfile(niigzsub,fullfile(sub3Ddir,[xsub{sub} '.nii.gz']));
        gunzip(fullfile(sub3Ddir,[xsub{sub} '.nii.gz']));
        delete(fullfile(sub3Ddir,[xsub{sub} '.nii.gz']));
        conversion = 0;
    elseif  numel(niftistart)==1
        conversion = 0;
    end
    
    niftistart = dir(niftisub);
    if conversion == -1 && numel(niftistart)~=1
    % Try DICOM (PARREC) conversion
        fils = dir(fullfile(subdcmdir,'*')); fils = fils(~[fils.isdir]);
        if numel(fils)>100
            conversion = 1;
        else
            fils = dir(fullfile(subdcmdir,'*.*')); fils = fils(~[fils.isdir]);
        end
        if numel(fils)>100
            conversion = 1;
        elseif numel(dir(fullfile(subdcmdir,'*.PAR')))==1 ||...
                numel(dir(fullfile(subdcmdir,'*.par')))==1
            conversion = 1;
        else
            fprintf('No viable images found for subject %s. Skipping...\n', xsub{sub});
        end
    end
    
    if conversion >= 0
        try
            se_subPreprocessing3D(subjects{sub},xsub{sub},structurals{1},conversion)

            isda = numel(dir(fullfile(sub3Ddir,'sm0wrp*.*'))) ...
                + numel(dir(fullfile(sub3Ddir,'smwrp*.*'))) ...
                + numel(dir(fullfile(sub3Ddir,'m0wrp*.*'))) ...
                + numel(dir(fullfile(sub3Ddir,'mwrp*.*')));
            if isda ~= 8
                fprintf (1,'%s\n',[xsub{sub} ': not (correctly) processed']);
            else
                fprintf (1,'%s\n',[xsub{sub} ': completed']);
            end

        catch err
            fprintf('Error processing %s: %s', xsub{sub}, err.message);

 	    if numel(fils) < 176 && numel(fils) >0
 	        fprintf(1,'%s\n',[xsub{sub} ': wrong number of DICOMS'])
 	    elseif numel(fils)==0
 	        fprintf(1,'%s\n',[xsub{sub} ': no 3D'])
 	    else
 	        fprintf(1,'%s\n',[xsub{sub} ' ... ups'])
 	    end
        end
    end
    
end

end

%%% ------ RESTING-STATE PREPROCESSING ------ %%%%%



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

fprintf('Found %i unprocessed of %i subjects with resting-state data..\n', numel(idx), numel(subjects));

if numel(idx) >= 1

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

end

