 clear, clc
 warning('off','all')

%--------------- readout of path-to-data ----------------------------
fid = fopen('_path2data','r'); datadir = fgetl(fid); fclose(fid);

Di = spm_select(Inf,'dir','Select subject directories','',fullfile(datadir,'DATA'));
 
% Temp directory for log and shell scripts
user_name = char(java.lang.System.getProperty('user.name')); % system('whoami');
user_name = strtrim(user_name);
tmp_dir = fullfile(pwd,'temp', 'dwi', user_name);
if exist(tmp_dir, 'dir'), rmdir(tmp_dir, 's'); end
mkdir(tmp_dir);


 remove_temp = true;
 
for i=1:size(Di,1)
    subjects{i} = strrep(Di(i,:),' ','');
    nx = strrep(Di(i,:),' ','');
    if nx(end) == filesep; nx(end) = []; end
    us = strfind(nx,filesep);
    nx = nx(us(end)+1:end);
    xsub{i} = nx;
end

[log_file,path] = uiputfile('*.*','Select log file');
log_file = fullfile(path,log_file);

if ~strcmp(log_file(end-3:end),'.log')
   log_file = [log_file '.log']; 
end

if exist(log_file, 'file')
    delete(log_file);
end

fix = fopen(log_file,'a+');
mkdir(fullfile(tmp_dir,'Shell2'));

fprintf('Subjects found: %i.\nStarting conversion...\n', numel(subjects));

for s=1:numel(subjects)
    
    if numel(strfind(computer,'MAC'))>0
        setenv('FSLDIR','/usr/local/fsl');
        setenv('FSLOUTPUTTYPE', 'NIFTI_GZ');
    end
    
    dwid = {'DTI','DTI60','DWI'};
    
    dcmdir = fullfile(subjects{s},'DICOM');
    dwidir = fullfile(subjects{s},'DWI');
    
    i = 1; fils = [];
    try
        while numel(fils)<60
            fils = dir(fullfile(dcmdir,dwid{i},'*'));
            if numel(fils)<60
                fils = dir(fullfile(dcmdir,dwid{i},'*.*'));
            end
            fils = fils(~[fils.isdir]); dwdir = dwid{i};
            i = i+1;
        end
    end
    
    conversion = -1;
    if numel(fils) >=60
       conversion = 1; 
    end
    % If we already have a NIFTI source with bvals & bvecs, copy this
    % instead of converting DICOMS (whether they exist or not)
    if numel(dir(fullfile(dcmdir,'DWI','*.nii.gz'))) > 0 && numel(dir(fullfile(dcmdir,'DWI','*.bval'))) > 0 && ...
            numel(dir(fullfile(dcmdir,'DWI','*.bvec'))) > 0
       conversion = 0;
    else %PARREC?
        
    end
     
    % Wenn DWI DICOMs da sind    
    %if numel(fils) >= 60
    if conversion > -1
        
        % Und noch nicht alles fertig ist
        if any([exist(fullfile(dwidir,'reg3G','Mean3G_warp2FA.nii.gz')) ...
                exist(fullfile(dwidir,[xsub{s} '_FA.nii.gz'])) ...
                exist(fullfile(dwidir,[xsub{s} '.bedpostX'],'mean_f1samples.nii.gz'))]~=2)
            
            shell_file = fullfile(tmp_dir,'Shell2',[xsub{s} '.sh']);
            
            % Alles mal weg, was da w?re
            [~, ~, ~] = rmdir(dwidir,'s');
            
            % DICOM conversion
%             mkdir(fullfile(dwidir))
%             mkdir(fullfile(dwidir,'DTI'))
            
            if conversion == 1
                % DICOM conversion
                mkdir(dwidir)
                mkdir(fullfile(dwidir,'DTI'))
                try
                    if numel(strfind(computer,'MAC'))>0
                        a = system([fullfile(pwd,'matfiles','dcm2niiMAC') ' -b ' fullfile(dcmdir,dwdir) ' ' fullfile(dcmdir,dwdir)]);
                    else
                        a = system([fullfile(pwd,'matfiles','dcm2nii') ' -b ' fullfile(dcmdir,dwdir) ' ' fullfile(dcmdir,dwdir)]);
                    end
                    fil = dir(fullfile(dcmdir,dwdir,'*.nii.gz'));
                    try
                        % Move converted files
                        movefile(fullfile(dcmdir,dwdir,fil(1).name),fullfile(dwidir,'raw.nii.gz')); 
                                    fil = dir(fullfile(dcmdir,dwdir,'*.bval'));
                        movefile(fullfile(dcmdir,dwdir,fil(1).name),fullfile(dwidir,'bvals')); 
                                    fil = dir(fullfile(dcmdir,dwdir,'*.bvec'));
                        movefile(fullfile(dcmdir,dwdir,fil(1).name),fullfile(dwidir,'bvecs')); 

                            % Remove any remaining nii, bvec, bval files
                        delete(fullfile(dcmdir,dwdir,'*.nii.gz'));
                        delete(fullfile(dcmdir,dwdir,'*.bval'));
                        delete(fullfile(dcmdir,dwdir,'*.bvec'));
                     catch err
                        fprintf(fix, err.message);
                     end
                catch

                    mkdir(fullfile(dwidir,'DICOM'))
                    copyfile(fullfile(dcmdir,dwdir),fullfile(dwidir,'DICOM'))
                    if numel(strfind(computer,'MAC'))>0
                        a = system([fullfile(pwd,'matfiles','dcm2niiMAC') ' -b ' fullfile(dwidir,'DICOM') ' ' fullfile(dwidir,'DICOM')]);
                    else
                        a = system([fullfile(pwd,'matfiles','dcm2nii') ' -b ' fullfile(dwidir,'DICOM') ' ' fullfile(dwidir,'DICOM')]);
                    end

                    fil = dir(fullfile(dwidir,'DICOM','*.nii.gz'));
                    movefile(fullfile(dwidir,'DICOM',fil(1).name),fullfile(dwidir,'raw.nii.gz'));
                    fil = dir(fullfile(dwidir,'DICOM','*.bval'));
                    movefile(fullfile(dwidir,'DICOM',fil(1).name),fullfile(dwidir,'bvals'));
                    fil = dir(fullfile(dwidir,'DICOM','*.bvec'));
                    movefile(fullfile(dwidir,'DICOM',fil(1).name),fullfile(dwidir,'bvecs'));

                    % Remove any remaining nii, bvec, bval files
                    delete(fullfile(dcmdir,dwdir,'*.nii.gz'));
                    delete(fullfile(dcmdir,dwdir,'*.bval'));
                    delete(fullfile(dcmdir,dwdir,'*.bvec'));

                    rmdir(fullfile(dwidir,'DICOM'),'s')
                end

            elseif conversion == 0
                mkdir(dwidir)
                mkdir(fullfile(dwidir,'DTI'))
                % Copy existing NIFTI + bvec + bval
                % Files must be named as subject ID
                copyfile(fullfile(dcmdir,'DWI',[xsub{s} '.nii.gz']),fullfile(dwidir,'raw.nii.gz'));
                copyfile(fullfile(dcmdir,'DWI',[xsub{s} '.bval']),fullfile(dwidir,'bvals'));
                copyfile(fullfile(dcmdir,'DWI',[xsub{s} '.bvec']),fullfile(dwidir,'bvecs'));
                
            end
            
            bvals = load(fullfile(dwidir,'bvals'));
            Q = find(bvals<100)-1;
            
            fid = fopen(shell_file,'w');
            fprintf(fid,'%s\n','#!/bin/bash');
            
            if numel(strfind(computer,'MAC'))>0
                fprintf(fid,'%s\n','FSLDIR=/usr/local/fsl');
                fprintf(fid,'%s\n','PATH=${FSLDIR}/bin:${PATH}');
                fprintf(fid,'%s\n','. ${FSLDIR}/etc/fslconf/fsl.sh');
                fprintf(fid,'%s\n','export FSLDIR PATH');
            end
            
            
            % Basic stuff until DTI
            fprintf(fid,'%s\n',[getenv('FSLDIR') '/bin/eddy_correct ' fullfile(dwidir,'raw.nii.gz') ' ' fullfile(dwidir,'data.nii.gz') ' 0']);
            
            for i=1:numel(Q)
                fprintf(fid,'%s\n',[getenv('FSLDIR') '/bin/fslroi ' fullfile(dwidir,'data.nii.gz') ' ' ...
                    fullfile(dwidir,['xno' int2str(i) '.nii.gz']) ' ' int2str(Q(i)) ' 1']);
            end
            fprintf(fid,'%s\n',[getenv('FSLDIR') '/bin/fslmerge -t ' fullfile(dwidir,['nodif.nii.gz']) ' ' fullfile(dwidir,['xno*.nii.gz'])]);
            fprintf(fid,'%s\n',[getenv('FSLDIR') '/bin/fslmaths ' fullfile(dwidir,['nodif.nii.gz']) ' -Tmean ' fullfile(dwidir,['nodif.nii.gz'])]);
            fprintf(fid,'%s\n',[getenv('FSLDIR') '/bin/bet ' fullfile(dwidir,['nodif.nii.gz']) ' '... 
                fullfile(dwidir,['nodif_brain.nii.gz']) ' -m -f 0.3 -R']);
            fprintf(fid,'%s\n',[getenv('FSLDIR') '/bin/dtifit -k ' fullfile(dwidir,'data.nii.gz') ...
                ' -m ' fullfile(dwidir,['nodif_brain_mask.nii.gz']) ...
                ' -r ' fullfile(dwidir,'bvecs') ...
                ' -b ' fullfile(dwidir,'bvals') ...
                ' -o ' fullfile(dwidir,'DTI','dti')]);
            
            fprintf(fid,'%s\n',['rm ' fullfile(dwidir,['nodif_brain.nii.gz'])]);
            fprintf(fid,'%s\n',['rm ' fullfile(dwidir,['xno*.nii.gz'])]);
            fprintf(fid,'%s\n',['rm ' fullfile(dwidir,['data.ecclog'])]);
            fprintf(fid,'%s\n',['cp ' fullfile(dwidir,'DTI','dti_FA.nii.gz') ' ' ...
                fullfile(dwidir,[xsub{s} '_FA.nii.gz'])]);
            
            % BEDPOST
            fprintf(fid,'%s\n',['mkdir ' fullfile(dwidir,[xsub{s}])]);
            fprintf(fid,'%s\n',['mv ' fullfile(dwidir,'nodif_brain_mask.nii.gz') ' ' ...
                fullfile(dwidir,[xsub{s}],'nodif_brain_mask.nii.gz')]);
            fprintf(fid,'%s\n',['mv ' fullfile(dwidir,'bvals') ' ' fullfile(dwidir,[xsub{s}],'bvals')]);
            fprintf(fid,'%s\n',['mv ' fullfile(dwidir,'bvecs') ' ' fullfile(dwidir,[xsub{s}],'bvecs')]);
            fprintf(fid,'%s\n',['mv ' fullfile(dwidir,'data.nii.gz') ' ' fullfile(dwidir,[xsub{s}],'data.nii.gz')]);
            
            fprintf(fid,'%s\n',[getenv('FSLDIR') '/bin/bedpostx ' fullfile(dwidir,xsub{s})]);
            
            fprintf(fid,'%s\n',['mv ' fullfile(dwidir,[xsub{s}],'*.*') ' ' fullfile(dwidir,[xsub{s} '.bedpostX'])]);
            fprintf(fid,'%s\n',['mv ' fullfile(dwidir,[xsub{s}],'*') ' ' fullfile(dwidir,[xsub{s} '.bedpostX'])]);
            
            
            % Warp to Mean3G
            fprintf(fid,'%s\n',['mkdir ' fullfile(dwidir,'reg3G')]);
            fprintf(fid,'%s\n',[getenv('FSLDIR') '/bin/flirt -in ' fullfile(dwidir,'DTI','dti_FA.nii.gz') ...
                ' -ref '  fullfile(pwd,'MaskenEtc','forDWI','Mean3G.nii.gz') ...
                ' -out '  fullfile(dwidir,'reg3G','FA_lin2Mean3G.nii.gz') ...
                ' -omat ' fullfile(dwidir,'reg3G','FA_lin2Mean3G.mat') ...
                ' -bins 256 -cost corratio -searchrx -180 180 -searchry -180 180 -searchrz -180 180 -dof 12 -interp spline']);
            
            fprintf(fid,'%s\n',[getenv('FSLDIR') '/bin/fnirt --in=' fullfile(dwidir,'DTI','dti_FA.nii.gz') ...
                ' --ref='     fullfile(pwd,'MaskenEtc','forDWI','Mean3G.nii.gz') ...
                ' --refmask=' fullfile(pwd,'MaskenEtc','forDWI','Mask.nii') ...
                ' --aff='    fullfile(dwidir,'reg3G','FA_lin2Mean3G.mat') ...
                ' --config=' fullfile(pwd,'G3_FAnormalization.cnf') ...
                ' --cout='   fullfile(dwidir,'reg3G','FA_warp2Mean3G.nii.gz') ...
                ' --iout='   fullfile(dwidir,'reg3G','FA_nlin2Mean3G.nii.gz')]);
            
            fprintf(fid,'%s\n',[getenv('FSLDIR') '/bin/invwarp '  ...
                ' --warp=' fullfile(dwidir,'reg3G','FA_warp2Mean3G.nii.gz') ...
                ' --out='  fullfile(dwidir,'reg3G','Mean3G_warp2FA.nii.gz')  ...
                ' --ref='  fullfile(dwidir,'DTI','dti_FA.nii.gz')]);
            
            
            % cleanup
            fprintf(fid,'%s\n',['rm -rf ' fullfile(dwidir,xsub{s})]);
            
            fclose(fid);
            fileattrib(shell_file,'+x');
            
            fprintf(fix,'%s\n',[fullfile(subjects{s}) ': setup shell file -> ' fullfile(pwd,'Shell',[xsub{s} '.sh'])]);
            %toDo{end+1} = fullfile(pwd,'Shell',[xsub{s} '.sh']);
        else
            
            % Cleanup ?
            try; delete(fullfile(dwidir,xsub{s},'*')); end
            try; delete(fullfile(dwidir,xsub{s},'.DS_Store')); end
            try; delete(fullfile(dwidir,'bvals')); end
            try; delete(fullfile(dwidir,'bvecs')); end
            try; delete(fullfile(dwidir,['raw.nii.gz'])); end
            try; [~, ~] = rmdir(fullfile(dwidir,[xsub{s}])); end
            
            fprintf(fix,'%s\n',[fullfile(subjects{s}) ': already done']);
        end
    else
        fprintf(fix,'%s\n',[fullfile(subjects{s}) ': Could not find 60 DICOM images']);
    end
end

fprintf('Done conversions. Starting parallel processing of scripts.\n');


try
	fprintf(fix, '--- Starting parallel processing of scripts ---\n');
	parfor i = 1 : numel(subjects);
	    lid = fopen(fullfile(tmp_dir,['parallel_' xsub{i} '.log']),'a+');
	    fprintf(lid, 'Starting script for subject %s.\n', xsub{i});
	    path_i = fullfile(tmp_dir,'Shell2',[xsub{i} '.sh']);
	    system(path_i);
	    try 
		mkdir(fullfile(subjects{i}, 'DWI', 'scripts')); 
	    catch err 
		fprintf(lid, 'Error creating scripts directory at %s', fullfile(subjects{i}, 'DWI', 'scripts'));
	    end
	    try 
		movefile(path_i,fullfile(subjects{i}, 'DWI', 'scripts', [xsub{i} '.sh'])); 
	    catch
		fprintf(lid, 'Error moving script to %s', fullfile(subjects{i}, 'DWI', 'scripts', [xsub{i} '.sh']));
	    end
	    fprintf(lid, 'Done script for subject %s.\n', xsub{i});
	    fclose(lid);
        try movefile(fullfile(tmp_dir,['parallel_' xsub{i} '.log']), ...
                fullfile(subjects{i}, 'DWI', 'logs', ['parallel_' xsub{i} '.log'])); catch; end
	end
	fprintf(fix, '--- Done parallel processing of scripts ---\n');
	
catch err
	fprintf(fix, 'Error in processing of parallel processing of scripts:\n');
	fprintf(fix, err.message);
	fprintf(err.message);
end


fclose(fix);

% Remove shell dir
if remove_temp
    rmdir(tmp_dir,'s');
end

fprintf('Done parallel processing of scripts.\n');

% clc;