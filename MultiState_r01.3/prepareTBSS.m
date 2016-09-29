clear;

fid = fopen('_path2data','r'); datadir = fgetl(fid); fclose(fid); cwd=pwd;

gname = spm_select(1,'mat','Select lookup file','',fullfile(cwd,'Lookups'));
load(gname)
gname = spm_str_manip(gname,'rt')

if numel(dir(fullfile('/usr/local/fsl','bin','fsl'))) == 1
    FSL = '/usr/local/fsl'; newtrax = 'probtrackx2_macosx'; % use on MAC
    setenv('FSLDIR',FSL); setenv('FSLOUTPUTTYPE', 'NIFTI_GZ');
elseif  numel(dir(fullfile('/usr/share/fsl/5.0','bin','fsl'))) == 1
    FSL = '/usr/share/fsl/5.0'; newtrax = 'probtrackx2_linux'; % use on cluster
    setenv('FSLDIR',FSL); setenv('FSLOUTPUTTYPE', 'NIFTI_GZ');
else
    error('Could not locate FSL')
end

mkdir(fullfile(cwd,'MaskenEtc','forDWI',gname,'NoDif'))
for i=1:numel(sub)
    if hasDTI(i)
        if     numel(dir(fullfile(datadir,'DATA',SubDir{i},sub{i},'DWI','reg3G','FA_nlin2Mean3G.nii.gz')))==1
            system([getenv('FSLDIR') '/bin/applywarp -i ' fullfile(datadir,'DATA',SubDir{i},sub{i},'DWI','nodif.nii.gz') ....
                ' -o ' fullfile(cwd,'MaskenEtc','forDWI',gname,'NoDif',sub{i}) ...
                ' -r ' fullfile(cwd,'MaskenEtc','forDWI','Mean3G.nii.gz') ...
                ' -w ' fullfile(datadir,'DATA',SubDir{i},sub{i},'DWI','reg3G','FA_warp2Mean3G.nii.gz')]); 
        elseif numel(dir(fullfile(datadir,'DATA',SubDir{i},sub{i},'DWI','reg3G','FA_nlin2Mean2G.nii.gz')))==1
            system([getenv('FSLDIR') '/bin/applywarp -i ' fullfile(datadir,'DATA',SubDir{i},sub{i},'DWI','nodif.nii.gz') ....
                ' -o ' fullfile(cwd,'MaskenEtc','forDWI',gname,'NoDif',sub{i}) ...
                ' -r ' fullfile(cwd,'MaskenEtc','forDWI','Mean3G.nii.gz') ...
                ' -w ' fullfile(datadir,'DATA',SubDir{i},sub{i},'DWI','reg3G','FA_warp2Mean2G.nii.gz')]); 
        else
            fprintf(' %i seems unprocessed.\n', sub{i});
        end
    end
end
system([getenv('FSLDIR') '/bin/fslmerge -t ' fullfile(cwd,'MaskenEtc','forDWI',gname,'all_NoDif.nii.gz') ' ' fullfile(cwd,'MaskenEtc','forDWI',gname,'NoDif','*.nii.gz') ]);
system([getenv('FSLDIR') '/bin/fslmaths ' fullfile(cwd,'MaskenEtc','forDWI',gname,'all_NoDif.nii.gz') ' -Tmean ' fullfile(cwd,'MaskenEtc','forDWI',gname,'mean_NoDif.nii.gz')  ]);




mkdir(fullfile(cwd,'MaskenEtc','forDWI',gname,'FA'))
mkdir(fullfile(cwd,'MaskenEtc','forDWI',gname,'Skeletonized'))
for i=1:numel(sub)
    if hasDTI(i)
        if numel(dir(fullfile(datadir,'DATA',SubDir{i},sub{i},'DWI','reg3G','FA_nlin2Mean3G.nii.gz')))==1
            copyfile(fullfile(datadir,'DATA',SubDir{i},sub{i},'DWI','reg3G','FA_nlin2Mean3G.nii.gz'),...
                fullfile(cwd,'MaskenEtc','forDWI',gname,'FA',[sub{i} '.nii.gz']))
        elseif numel(dir(fullfile(datadir,'DATA',SubDir{i},sub{i},'DWI','reg3G','FA_nlin2Mean2G.nii.gz')))==1
            copyfile(fullfile(datadir,'DATA',SubDir{i},sub{i},'DWI','reg3G','FA_nlin2Mean2G.nii.gz'),...
                fullfile(cwd,'MaskenEtc','forDWI',gname,'FA',[sub{i} '.nii.gz']))
        else
            fprintf(' %i needs preprocessing!?\n', sub{i});
        end
    end
end

system([getenv('FSLDIR') '/bin/fslmerge -t ' fullfile(cwd,'MaskenEtc','forDWI',gname,'all_FA.nii.gz') ' ' fullfile(cwd,'MaskenEtc','forDWI',gname,'FA','*.nii.gz') ]);
system([getenv('FSLDIR') '/bin/fslmaths ' fullfile(cwd,'MaskenEtc','forDWI',gname,'all_FA.nii.gz') ' -Tmean ' fullfile(cwd,'MaskenEtc','forDWI',gname,'meanFA.nii.gz')  ]);

system([getenv('FSLDIR') '/bin/bet ' fullfile(cwd,'MaskenEtc','forDWI',gname,'meanFA.nii.gz') ' ' fullfile(cwd,'MaskenEtc','forDWI',gname,'meanFA_brain.nii.gz')  '  -m -n -R']);
system([getenv('FSLDIR') '/bin/fslmaths ' fullfile(cwd,'MaskenEtc','forDWI',gname,'meanFA_brain_mask.nii.gz') ' -eroF -eroF -eroF ' fullfile(cwd,'MaskenEtc','forDWI',gname,'meanFA_brain_mask.nii.gz')  ]);
system([getenv('FSLDIR') '/bin/fslmaths ' fullfile(cwd,'MaskenEtc','forDWI',gname,'meanFA.nii.gz') ' -mas ' fullfile(cwd,'MaskenEtc','forDWI',gname,'meanFA_brain_mask.nii.gz') ' -thr .1  -bin '  fullfile(cwd,'MaskenEtc','forDWI',gname,'mean_FA_mask.nii.gz') '  -odt char']);
system([getenv('FSLDIR') '/bin/fslmaths ' fullfile(cwd,'MaskenEtc','forDWI',gname,'all_FA.nii.gz') ' -mas ' fullfile(cwd,'MaskenEtc','forDWI',gname,'mean_FA_mask.nii.gz') ' '  fullfile(cwd,'MaskenEtc','forDWI',gname,'all_FA.nii.gz')]);
system([getenv('FSLDIR') '/bin/fslmaths ' fullfile(cwd,'MaskenEtc','forDWI',gname,'all_FA.nii.gz') ' -Tmean ' fullfile(cwd,'MaskenEtc','forDWI',gname,'mean_FA.nii.gz')]);

system([getenv('FSLDIR') '/bin/tbss_skeleton -i ' fullfile(cwd,'MaskenEtc','forDWI',gname,'mean_FA.nii.gz') ' -o ' fullfile(cwd,'MaskenEtc','forDWI',gname,'mean_FA_skeleton.nii.gz')]);

system([getenv('FSLDIR') '/bin/fslmaths ' fullfile(cwd,'MaskenEtc','forDWI',gname,'mean_FA_skeleton.nii.gz') ' -thr .2 -bin  ' fullfile(cwd,'MaskenEtc','forDWI',gname,'mean_FA_skeleton_mask.nii.gz')]);
system([getenv('FSLDIR') '/bin/fslmaths ' fullfile(cwd,'MaskenEtc','forDWI',gname,'mean_FA_mask.nii.gz') ' -mul -1 -add 1 -add  ' fullfile(cwd,'MaskenEtc','forDWI',gname,'mean_FA_skeleton_mask.nii.gz') ' ' fullfile(cwd,'MaskenEtc','forDWI',gname,'mean_FA_skeleton_mask_dst.nii.gz')]);
system([getenv('FSLDIR') '/bin/distancemap -i ' fullfile(cwd,'MaskenEtc','forDWI',gname,'mean_FA_skeleton_mask_dst.nii.gz') ' -o ' fullfile(cwd,'MaskenEtc','forDWI',gname,'mean_FA_skeleton_mask_dst.nii.gz')]);

for i=1:numel(sub)
    if numel(dir(fullfile('/usr/local/fsl','bin','fsl'))) == 1
        FSL = '/usr/local/fsl'; newtrax = 'probtrackx2_macosx'; % use on MAC
        setenv('FSLDIR',FSL); setenv('FSLOUTPUTTYPE', 'NIFTI_GZ');
    elseif  numel(dir(fullfile('/usr/share/fsl/5.0','bin','fsl'))) == 1
        FSL = '/usr/share/fsl/5.0'; newtrax = 'probtrackx2_linux'; % use on cluster
        setenv('FSLDIR',FSL); setenv('FSLOUTPUTTYPE', 'NIFTI_GZ');
    else
        error('Could not locate FSL')
    end
    system([getenv('FSLDIR') '/bin/tbss_skeleton -i ' fullfile(cwd,'MaskenEtc','forDWI',gname,'mean_FA.nii.gz') ' -p .2 ' fullfile(cwd,'MaskenEtc','forDWI',gname,'mean_FA_skeleton_mask_dst.nii.gz') ' ' ...
        getenv('FSLDIR') '/data/standard/LowerCingulum_1mm ' fullfile(cwd,'MaskenEtc','forDWI',gname,'FA',[sub{i} '.nii.gz']) ' ' ...
        fullfile(cwd,'MaskenEtc','forDWI',gname,'Skeletonized',[sub{i} '.nii.gz'])]);
end

rmdir(fullfile(cwd,'MaskenEtc','forDWI',gname,'FA'),'s')
delete(fullfile(cwd,'MaskenEtc','forDWI',gname,'meanFA_brain.nii.gz'))
delete(fullfile(cwd,'MaskenEtc','forDWI',gname,'meanFA_brain_mask.nii.gz'))
delete(fullfile(cwd,'MaskenEtc','forDWI',gname,'all_FA.nii.gz'))
delete(fullfile(cwd,'MaskenEtc','forDWI',gname,'mean_FA_mask.nii.gz'))
delete(fullfile(cwd,'MaskenEtc','forDWI',gname,'mean_FA_skeleton.nii.gz'))

