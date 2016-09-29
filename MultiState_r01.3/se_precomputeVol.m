function se_precomputeVol(des)

global NetFlags


fid = fopen('_path2data','r'); datadir = fgetl(fid); fclose(fid);

[~, ~, ~] = mkdir(fullfile(pwd,'Projects','VBM','Regional',spm_str_manip(NetFlags.Group,'rt')));
fid = fopen(fullfile(pwd,'Projects','VBM','Regional',spm_str_manip(NetFlags.Group,'rt'),[NetFlags.project '_' des '_subjects.txt']),'w+');
fprintf(fid,'%s\t%s\t','Subject','Patient');
for i=1:numel(NetFlags.Covariates)
    fprintf(fid,'%s\t',NetFlags.Covariates{i});
end
fprintf(fid,'%s\t','is Patient');
for ii=1:numel(NetFlags.sub)
    fprintf(fid,'\n%s\t',NetFlags.sub{ii});
    fprintf(fid,'%4.0f\t',NetFlags.isPat(ii));
    
    for i=1:numel(NetFlags.Covariates)
        fprintf(fid,'%4.1f\t',NetFlags.Cov(ii,i));
    end
end
fclose(fid);


NN = [];

z = nan(numel(NetFlags.sub),numel(NetFlags.Label));

clear VOI


if NetFlags.Coordinates == 0
    
    if get(NetFlags.CN,'Value')==2
        try
            for ana = 1:numel(NetFlags.Label)
                VOI(ana) = spm_vol(fullfile(pwd,'Projects','VBM','Regional',spm_str_manip(NetFlags.Group,'rt'),NetFlags.project,['w' spm_str_manip(NetFlags.VOI{ana},'t')]));
            end
        catch
            [~, ~, ~] = mkdir(fullfile(pwd,'Projects','VBM','Regional',spm_str_manip(NetFlags.Group,'rt'),NetFlags.project));
            load(fullfile(pwd,'matfiles','VOIdeformations.mat'));
            matlabbatch{1}.spm.util.defs.comp{1}.sn2def.matname{1} = fullfile(pwd,'MaskenEtc','forVBM',['Colin2' spm_str_manip(NetFlags.Group,'rt') '_sn.mat']);
            matlabbatch{1}.spm.util.defs.savedir.saveusr{1} = fullfile(pwd,'Projects','VBM','Regional',spm_str_manip(NetFlags.Group,'rt'),NetFlags.project);
            matlabbatch{1}.spm.util.defs.fnames = NetFlags.VOI;
            spm_jobman('run',matlabbatch)
            
            clear VOI
            for ana = 1:numel(NetFlags.Label)
                VOI(ana) = spm_vol(fullfile(pwd,'Projects','VBM','Regional',spm_str_manip(NetFlags.Group,'rt'),NetFlags.project,['w' spm_str_manip(NetFlags.VOI{ana},'t')]));
            end
            
        end
    elseif get(NetFlags.CN,'Value')==3
        try
            for ana = 1:numel(NetFlags.Label)
                VOI(ana) = spm_vol(fullfile(pwd,'Projects','VBM','Regional',spm_str_manip(NetFlags.Group,'rt'),NetFlags.project,['w' spm_str_manip(NetFlags.VOI{ana},'t')]));
            end
        catch
            [~, ~, ~] = mkdir(fullfile(pwd,'Projects','VBM','Regional',spm_str_manip(NetFlags.Group,'rt'),NetFlags.project));
            load(fullfile(pwd,'matfiles','VOIdeformations.mat'));
            matlabbatch{1}.spm.util.defs.comp{1}.sn2def.matname{1} = fullfile(pwd,'MaskenEtc','forVBM',['MNI152T12' spm_str_manip(NetFlags.Group,'rt') '_sn.mat']);
            matlabbatch{1}.spm.util.defs.savedir.saveusr{1} = fullfile(pwd,'Projects','VBM','Regional',spm_str_manip(NetFlags.Group,'rt'),NetFlags.project);
            matlabbatch{1}.spm.util.defs.fnames = NetFlags.VOI;
            spm_jobman('run',matlabbatch)
            
            clear VOI
            for ana = 1:numel(NetFlags.Label)
                VOI(ana) = spm_vol(fullfile(pwd,'Projects','VBM','Regional',spm_str_manip(NetFlags.Group,'rt'),NetFlags.project,['w' spm_str_manip(NetFlags.VOI{ana},'t')]));
            end
            
        end
        
        
    else
        for ana = 1:numel(NetFlags.Label)
            try
                VOI(ana) = spm_vol(NetFlags.VOI{ana});
            catch
                VOI(ana) = spm_vol(fullfile(pwd,'VOIs','VOIfiles',spm_str_manip(NetFlags.VOI{ana},'t')));
            end
        end
        
    end
    
    
    for ana = 1:numel(NetFlags.Label)
        dat = spm_read_vols(VOI(ana));
        tmp = find(dat); Zdat{ana} = dat(tmp); clear XYZ
        [XYZ(:,1) XYZ(:,2) XYZ(:,3)] = ind2sub(VOI(ana).dim,tmp);
        tXYZ{ana} = (VOI(ana).mat) * [XYZ'; ones(1,size(XYZ,1))];
        dete(ana) = abs(det(VOI(ana).mat));
    end
    
else
Vgm = spm_vol(fullfile(pwd,'MaskenEtc','forVBM',[spm_str_manip(NetFlags.Group,'rt') '_c1meanT1.nii']));
gm  = spm_read_vols(Vgm);
ind = find(gm>.1); clear XYZ
[XYZ(:,1) XYZ(:,2) XYZ(:,3)] = ind2sub(Vgm.dim,ind); clear ind gm
XYZ = (Vgm.mat) * [XYZ'; ones(1,size(XYZ,1))];
    for ana = 1:numel(NetFlags.Label)
        D = sqrt((XYZ(1,:)-NetFlags.X(ana)).^2+(XYZ(2,:)-NetFlags.Y(ana)).^2+(XYZ(3,:)-NetFlags.Z(ana)).^2);
        tXYZ{ana} = XYZ(:,D<=5);
        dete(ana) = 1;
        Zdat{ana} = normpdf(D(D<=5),0,1.5); Zdat{ana} = Zdat{ana}/sum(Zdat{ana}); Zdat{ana} = Zdat{ana}';
    end
end

pref = {'m','m0','sm','sm0'};

for s=1:numel(NetFlags.sub)
    fprintf(1,'%s\n',['Subject ' int2str(s) ' / ' int2str(numel(NetFlags.sub))]);
    
    
    fil = dir(fullfile(datadir,'DATA',NetFlags.SubDir{s},NetFlags.sub{s},'3D',[pref{get(NetFlags.mtype,'value')} 'wrp1*.nii']));
    V1 = spm_vol(fullfile(datadir,'DATA',NetFlags.SubDir{s},NetFlags.sub{s},'3D',fil(1).name));
    try
        fil = dir(fullfile(datadir,'DATA',NetFlags.SubDir{s},NetFlags.sub{s},'3D',[pref{get(NetFlags.mtype,'value')} 'wrp2*.nii']));
        V2 = spm_vol(fullfile(datadir,'DATA',NetFlags.SubDir{s},NetFlags.sub{s},'3D',fil(1).name));
    end
    
    if get(NetFlags.Mask,'value')==1
        V = V1;
    else
        V = V2;
    end
    clear vol
    
    for ana = 1:numel(NetFlags.Label)
        vol(ana) = nansum(spm_get_data(V,inv(V.mat) * tXYZ{ana}).*Zdat{ana}')*dete(ana);
    end
    %     if get(NetFlags.Glob,'value')==1
    %         tv = spm_read_vols(V); tv = nansum(tv(:)); tv = tv*abs(det(V.mat));
    %         z(s,:) = vol/tv;
    %     elseif get(NetFlags.Glob,'value')==2
    %         tv = spm_read_vols(V1)+spm_read_vols(V2); tv = nansum(tv(:)); tv = tv*abs(det(V.mat));
    %         z(s,:) = vol/tv;
    %     else
    z(s,:) = vol;
    %    end
end




NetFlags.z = z;

Label = NetFlags.Label';
save(fullfile(pwd,'Projects','VBM','Regional',spm_str_manip(NetFlags.Group,'rt'),[NetFlags.project '_' des '.mat']),'z','Label');





wrz = [];
wn  = {};
cnt = 1;

for i=1:numel(Label)
    wrz(cnt,:) = z(:,i)';
    wn{cnt} = [Label{i}];
    cnt = cnt+1;
end
wrz = wrz';

fid = fopen(fullfile(pwd,'Projects','VBM','Regional',spm_str_manip(NetFlags.Group,'rt'),[NetFlags.project '_' des '_volumes.txt']),'w+');
fprintf(fid,'%s\t','Subject');
for i=1:numel(wn)
    fprintf(fid,'%s\t',wn{i});
end
for ii=1:numel(NetFlags.sub)
    fprintf(fid,'\n%s\t',NetFlags.sub{ii});
    
    for i=1:numel(wn)
        fprintf(fid,'%10.5f\t',wrz(ii,i));
    end
end
fclose(fid);



