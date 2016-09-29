function se_computeTBSS

global GLMFlags

mgroup = spm_str_manip(GLMFlags.Group,'rt');
warCov = size(GLMFlags.rawCov,2);

fid = fopen('_path2data','r'); datadir = fgetl(fid); fclose(fid);


analyses  = {'GroupDifferences','GroupDifferences_Covariate','Group_x_Covariate','All_Covariates','Patients_Covariates'};

if isfield(GLMFlags,'Label') && numel(GLMFlags.Label)>0 && any(GLMFlags.ActiveCov>warCov)
    folderNow = spm_input('Title','+1','s',GLMFlags.project);
else
    folderNow = spm_input('Title','+1','s',[analyses{get(GLMFlags.Analysis,'value')}]);

end


datafolder = [];
if get(GLMFlags.Confx,'Value')>1
    
    Conf      = [get(GLMFlags.Conf1,'Value') get(GLMFlags.Conf2,'Value') get(GLMFlags.Conf3,'Value')]; Conf = unique(Conf(Conf>1))-1; cont = [];
    for i=1:numel(Conf);
        datafolder = [datafolder GLMFlags.Covariates{Conf(i)}(1:2) '_' ];
    end
    datafolder = [datafolder  int2str(get(GLMFlags.Confx,'Value')-1) '_'];
else
    datafolder = 'NoConf_';
end

folderNow = [datafolder '_' folderNow];
datafolder = [datafolder 'DATA'];


[~,~,~] = mkdir(fullfile(pwd,'Projects','TBSS','GLM',mgroup,datafolder));
[~,~,~] = mkdir(fullfile(pwd,'Projects','TBSS','GLM',mgroup,folderNow));




OK = 1;
for s=1:numel(GLMFlags.sub)
    if ~exist(fullfile(pwd,'Projects','TBSS','GLM',mgroup,datafolder,[GLMFlags.sub{s} '.nii']));
        OK=0;
    end
end


gunzip(fullfile(pwd,'MaskenEtc','forDWI',spm_str_manip(GLMFlags.Group,'rt'),'mean_FA_skeleton_mask.nii.gz'))
Vgm = spm_vol(fullfile(pwd,'MaskenEtc','forDWI',spm_str_manip(GLMFlags.Group,'rt'),'mean_FA_skeleton_mask.nii'));
gm  = spm_read_vols(Vgm);
ind = find(gm>.01); clear XYZ
[XYZ(:,1) XYZ(:,2) XYZ(:,3)] = ind2sub(Vgm.dim,ind); clear ind gm





if OK == 0
    
    dat = zeros(numel(GLMFlags.sub),size(XYZ,1));
    for s=1:numel(GLMFlags.sub)
        gunzip(fullfile(pwd,'MaskenEtc','forDWI',spm_str_manip(GLMFlags.Group,'rt'),'Skeletonized',[GLMFlags.sub{s} '.nii.gz']))
        V        = spm_vol(fullfile(pwd,'MaskenEtc','forDWI',spm_str_manip(GLMFlags.Group,'rt'),'Skeletonized',[GLMFlags.sub{s} '.nii']));
        sXYZ     = V.mat \ Vgm.mat * [XYZ'; ones(1,size(XYZ,1))];
        dat(s,:) = spm_sample_vol(V,XYZ(:,1),XYZ(:,2),XYZ(:,3),1);
        delete(fullfile(pwd,'MaskenEtc','forDWI',spm_str_manip(GLMFlags.Group,'rt'),'Skeletonized',[GLMFlags.sub{s} '.nii']))
    end
    clear sXYZ
    
    
    
    
    if get(GLMFlags.Confx,'Value')>1
        
        c   = mean(dat);
        Conf = [get(GLMFlags.Conf1,'Value') get(GLMFlags.Conf2,'Value') get(GLMFlags.Conf3,'Value')]; Conf = unique(Conf(Conf>1))-1; cont = [];
        for iii=1:numel(Conf)
            if (any(round(GLMFlags.Cov(:,Conf(iii)))~=GLMFlags.Cov(:,Conf(iii))) & numel(unique(GLMFlags.Cov(:,Conf(iii))))>2) | ...
                    numel(unique(GLMFlags.Cov(:,Conf(iii))))>6; cont = [cont iii];
            end
        end
        dat = dat-repmat(c,size(dat,1),1);
        mod = mat2cell(GLMFlags.Cov(:,Conf),size(GLMFlags.Cov,1),[ones(1,numel(Conf))]);
        ord = get(GLMFlags.Confx,'Value')-1;
        useP = 0; try; parfor i=1:10; i; end; useP=1; end
        
        if useP
            parfor i=1:size(dat,2)
                [p,table,stats] = anovan(dat(:,i),mod,'model',ord,'display','off','sstype',3,'continuous',cont);
                dat(:,i) = stats.resid+c(i);
            end
        else
            for i=1:size(dat,2)
                [p,table,stats] = anovan(dat(:,i),mod,'model',ord,'display','off','sstype',3,'continuous',cont);
                dat(:,i) = stats.resid+c(i);
            end
        end
        
    end
    
    
    % Jetzt wieder rausschreiben
    for s=1:numel(GLMFlags.sub)
        Vo = rmfield(Vgm,'pinfo');
        Vo = rmfield(Vo,'private');
        Vo.fname = fullfile(pwd,'Projects','TBSS','GLM',mgroup,datafolder,[GLMFlags.sub{s} '.nii']);
        Vo        = spm_write_vol(Vo,accumarray(XYZ,dat(s,:),Vo.dim));
    end
    
end
clear dat





% Check if only Pat
%    matlabbatch{1}.spm.stats.factorial_design.des.mreg.scans(GLMFlags.isPat'==1 & sum(isnan(GLMFlags.Cov(:,GLMFlags.ActiveCov)),2)==0);

if get(GLMFlags.Confx,'Value')>1
    
    Conf = [get(GLMFlags.Conf1,'Value') get(GLMFlags.Conf2,'Value') get(GLMFlags.Conf3,'Value')]; Conf = unique(Conf(Conf>1))-1; cont = [];
    for iii=1:numel(Conf)
        if (any(round(GLMFlags.Cov(:,Conf(iii)))~=GLMFlags.Cov(:,Conf(iii))) & numel(unique(GLMFlags.Cov(:,Conf(iii))))>2) | ...
                numel(unique(GLMFlags.Cov(:,Conf(iii))))>6; cont = [cont iii];
        end
    end
    for i=1:size(GLMFlags.rawCov,2)
        if ~any(Conf==i)
        Q = ~isnan(GLMFlags.rawCov(:,i));
        c = mean(GLMFlags.rawCov(Q,i));
        [p,table,stats] = anovan(GLMFlags.rawCov(Q,i)-c,mat2cell(GLMFlags.Cov(Q,Conf),sum(Q),[ones(1,numel(Conf))]),...
            'model',get(GLMFlags.Confx,'Value')-1,'display','off','sstype',3,'continuous',cont);
        GLMFlags.Cov(Q,i) = stats.resid+c;
        else
            GLMFlags.Cov(:,i) = GLMFlags.rawCov(:,i);
        end
    end
end




if isfield(GLMFlags,'Label') && numel(GLMFlags.Label)>0 && any(GLMFlags.ActiveCov>warCov)
    
    VBM = spm_vol(fullfile(pwd,'MaskenEtc','forVBM',[spm_str_manip(GLMFlags.Group,'rt') '_c1meanT1.nii']));
    vbm  = spm_read_vols(VBM);
    ind = find(vbm>.1); clear vXYZ
    [vXYZ(:,1) vXYZ(:,2) vXYZ(:,3)] = ind2sub(VBM.dim,ind); clear ind vbm
    XYZmm = (VBM.mat) * [vXYZ'; ones(1,size(vXYZ,1))];
    vois = GLMFlags.ActiveCov(find(GLMFlags.ActiveCov>warCov))-warCov;
    
    for ana = 1:numel(vois)
        if GLMFlags.Coordinates == 0
            VOI(ana) = spm_vol(GLMFlags.VOI{vois(ana)});
            dat = find(spm_read_vols(VOI(ana))); clear tmp
            [tmp(:,1) tmp(:,2) tmp(:,3)] = ind2sub(VOI(ana).dim,dat);
            sXYZ{ana} = tmp';
        else
            VOI(ana) = Vgm;
            D = sqrt((XYZmm(1,:)-GLMFlags.X(vois(ana))).^2+(XYZmm(2,:)-GLMFlags.Y(vois(ana))).^2+(XYZmm(3,:)-GLMFlags.Z(vois(ana))).^2);
            sXYZ{ana} = vXYZ(D<=5,:)';
        end
    end
    
    dat = zeros(numel(GLMFlags.sub),numel(vois)); try; do = get(GLMFlags.VBM,'Value'); catch; do=2; end
    for s=1:numel(GLMFlags.sub)
        if do == 2
            Vo = spm_vol(fullfile(datadir,'DATA',GLMFlags.SubDir{s},GLMFlags.sub{s},'3D',['mwrp1' GLMFlags.sub{s} '.nii']));
            for ana = 1:numel(vois)
                nowXYZ     = inv(Vo.mat) * VOI(ana).mat * [sXYZ{ana}; ones(1,size(sXYZ{ana},2))];
                dat(s,ana) = nansum(spm_get_data(Vo,nowXYZ));
            end
        elseif do == 3
            Vo = spm_vol(fullfile(datadir,'DATA',GLMFlags.SubDir{s},GLMFlags.sub{s},'3D',['m0wrp1' GLMFlags.sub{s} '.nii']));
            for ana = 1:numel(vois)
                nowXYZ     = inv(Vo.mat) * VOI(ana).mat * [sXYZ{ana}; ones(1,size(sXYZ{ana},2))];
                dat(s,ana) = nansum(spm_get_data(Vo,nowXYZ));
            end
        else
            GLMFlags.ActiveCov(GLMFlags.ActiveCov>warCov) = [];
        end
    end
    
    if do>1
        GLMFlags.Cov(:,GLMFlags.ActiveCov(GLMFlags.ActiveCov>warCov))        = dat;
    end
    
end






load(fullfile(pwd,'matfiles','CovVBM.mat'));
for s=1:numel(GLMFlags.sub)
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.scans{s} = ...
        fullfile(pwd,'Projects','TBSS','GLM',mgroup,datafolder,[GLMFlags.sub{s} '.nii']);
end


matlabbatch{1}.spm.stats.factorial_design.dir{1} = fullfile(pwd,'Projects','TBSS','GLM',mgroup,folderNow);


mcov = matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(1);
ncov = matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(2);
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov = matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(1);



cnt = 1;
cov = nan(1,numel(GLMFlags.sub));
if get(GLMFlags.Analysis,'value')<4
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(cnt) = mcov;
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(cnt).c = GLMFlags.isPat;
    cov(cnt,:) = matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(cnt).c;
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(cnt).cname = 'Main Effect';
    cnt = cnt+1;
end

if get(GLMFlags.Analysis,'value')==2 | get(GLMFlags.Analysis,'value')==4
    for i=1:numel(GLMFlags.ActiveCov)
        matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(cnt) = mcov;
        matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(cnt).c = GLMFlags.Cov(:,GLMFlags.ActiveCov(i));
        cov(cnt,:) = matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(cnt).c;
        matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(cnt).cname = GLMFlags.Covariates{GLMFlags.ActiveCov(i)};
        cnt = cnt+1;
    end
end

if get(GLMFlags.Analysis,'value')==3
    for i=1:numel(GLMFlags.ActiveCov)
        matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(cnt) = ncov;
        matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(cnt).c = GLMFlags.Cov(:,GLMFlags.ActiveCov(i));
        matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(cnt).c(GLMFlags.isPat==0) = 0;
        matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(cnt).c(GLMFlags.isPat==1) = ...;
        matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(cnt).c(GLMFlags.isPat==1) - ...
        mean(matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(cnt).c(GLMFlags.isPat==1));
        cov(cnt,:) = matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(cnt).c;
        matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(cnt).cname = [GLMFlags.Covariates{GLMFlags.ActiveCov(i)} ' Pat'];
        cnt = cnt+1;

        matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(cnt) = ncov;
        matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(cnt).c = GLMFlags.Cov(:,GLMFlags.ActiveCov(i));
        matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(cnt).c(GLMFlags.isPat==1) = 0;
        matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(cnt).c(GLMFlags.isPat==0) = ...;
        matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(cnt).c(GLMFlags.isPat==0) - ...
        mean(matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(cnt).c(GLMFlags.isPat==0));
        cov(cnt,:) = matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(cnt).c;
        matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(cnt).cname = [GLMFlags.Covariates{GLMFlags.ActiveCov(i)} ' Con'];
        cnt = cnt+1;

    end
end

if get(GLMFlags.Analysis,'value')==5
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.scans = ...
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.scans(GLMFlags.isPat'==1 & sum(isnan(GLMFlags.Cov(:,GLMFlags.ActiveCov)),2)==0);
    for i=1:numel(GLMFlags.ActiveCov)
        matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(cnt) = mcov;
        matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(cnt).c = GLMFlags.Cov(GLMFlags.isPat'==1 & sum(isnan(GLMFlags.Cov(:,GLMFlags.ActiveCov)),2)==0,GLMFlags.ActiveCov(i));
        matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(cnt).cname = [GLMFlags.Covariates{GLMFlags.ActiveCov(i)}];
        cnt = cnt+1;

    end
end



% Wie mit NaNs umgehen -> diese Raus, diese 0
%for i=1:cnt-1
%    matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(cnt).c
%end


cd(matlabbatch{1}.spm.stats.factorial_design.dir{1})
save Design matlabbatch
spm_jobman('run',matlabbatch)



