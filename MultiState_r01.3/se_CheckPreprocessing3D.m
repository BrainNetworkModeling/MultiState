clear, clc; warning off all

force = 0;

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

structurals     = {'3D'};

% First determine which subjects need processing, to ensure an even
% distribution across workers

idx={};
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
        	if numel(dir(fullfile(subjects{sub},structurals{1}))) > 0
	            idx(i) = sub;
	            i = i + 1;
	        end
        end

    end
end

if numel(idx) == 0
    fprintf('No subjects to process! Aborting.');
    return; 
end

fprintf('Found %i unprocessed of %i subjects..\n', numel(idx), numel(subjects));

subjects = subjects(idx);
xsub = xsub(idx);


parfor sub = 1:numel(subjects)
    
    % What conversion is necessary? 0 = none (NIFTI exists), 1 = DICOM
    conversion = -1; fils = [];

	subdcmdir = fullfile(subjects{sub},'DICOM',structurals{1});
	sub3Ddir  = fullfile(subjects{sub},structurals{1});
	niigzsub = fullfile(subdcmdir,[xsub{sub} '.nii.gz']);
	niftisub = fullfile(sub3Ddir,[xsub{sub} '.nii']);

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
            try movefile(fullfile(subdcmdir,niigznames(1).name),niigzsub); end
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
            se_subPreprocessing3D(subjects{sub},xsub{sub},structurals{1},conversion,force)

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
