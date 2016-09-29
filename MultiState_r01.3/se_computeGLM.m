function se_computeGLM

global GLMFlags

cwd=pwd;

fid = fopen('_path2data','r'); datadir = fgetl(fid); fclose(fid);

%%% ----    
try; GLMFlags.gotMask = get(GLMFlags.Mask,'value'); end
try; GLMFlags.gotPCA = get(GLMFlags.PCA,'value'); end
try; GLMFlags.gotGlob = get(GLMFlags.Glob,'value'); end
try; GLMFlags.gotFilt = get(GLMFlags.Filt,'value'); end
try; GLMFlags.gotBlock = get(GLMFlags.Block,'value'); end
try; GLMFlags.gotAnalysis = get(GLMFlags.Analysis,'value'); end
try; GLMFlags.gotVBM = get(GLMFlags.VBM,'Value'); end
try; GLMFlags.keep; catch; GLMFlags.keep = 0; end
try; GLMFlags.gotConf1 = get(GLMFlags.Conf1,'Value'); end  
try; GLMFlags.gotConf2 = get(GLMFlags.Conf2,'Value'); end
try; GLMFlags.gotConf3 = get(GLMFlags.Conf3,'Value'); end
try; GLMFlags.gotConf1t = get(GLMFlags.Conf1t,'Value'); end  
try; GLMFlags.gotConf2t = get(GLMFlags.Conf2t,'Value'); end
try; GLMFlags.gotConf3t = get(GLMFlags.Conf3t,'Value'); end
try; GLMFlags.gotConfx = get(GLMFlags.Confx,'Value'); end

GM   = {'NoGMmask','IndGMmask','GroupGMmask'};
PCA  = {'PCA','NoPCA','FIX'};
Glob = {'TSR','lTSR','WMCSF','lWMCSF','GSR','lGSR','NoGSR'};
FILT = {'BP','HP'};

if GLMFlags.Coordinates==0
    des = [GM{GLMFlags.gotMask} '_' PCA{GLMFlags.gotPCA} '_' Glob{GLMFlags.gotGlob} '_' FILT{GLMFlags.gotFilt}];
else
    des = [GM{GLMFlags.gotMask} '_' PCA{GLMFlags.gotPCA}  '_' Glob{GLMFlags.gotGlob} '_' FILT{GLMFlags.gotFilt}...
        '_' strrep(num2str(GLMFlags.Radius,'%2.1f'),'.','')];
end
%%% ----

if GLMFlags.gotPCA==3
    RSdir='RS_fix';
    Qreg = [];       % 24 movement regressors already corrected for in the FIX denoising
else
    RSdir='RS';
    Qreg = 9:32;    % 24 movement regressors from the realignment parameter
end

%        Nuisance regressors as estimated in the preprocessing:
%        Nr      4     4      6   6    6    6         5
%        No      1-    5-     9-  15-  21-  27-       33-    37
%        reg = [gx2' gx2'.^2 rp rp.^2 rpp rpp.^2 COEFF(:,1:5)];

if GLMFlags.gotPCA==1
    Qreg = [Qreg 33:37];
end
if GLMFlags.gotGlob==1
    Qreg = [Qreg 1 2 3 5 6 7];
elseif GLMFlags.gotGlob==2
    Qreg = [Qreg 1:3];
elseif GLMFlags.gotGlob==3
    Qreg = [Qreg 2 3 6 7];
elseif GLMFlags.gotGlob==4
    Qreg = [Qreg 2:3];
elseif GLMFlags.gotGlob==5
    Qreg = [Qreg 4 8];
elseif GLMFlags.gotGlob==6
    Qreg = [Qreg 4];
end

if GLMFlags.Filt==1
    Filt = [0.01, 0.08];
else
    Filt = [0.01, inf];
end    
      

Conf = [GLMFlags.gotConf1 GLMFlags.gotConf2 GLMFlags.gotConf3];
Conft = [GLMFlags.gotConf1t GLMFlags.gotConf2t GLMFlags.gotConf3t];

[Conf ia] = unique(Conf(Conf>1));
Conf = Conf-1; cont = find((Conft(ia)-1)>0);

modl = GLMFlags.gotConfx-1;

if modl>0 & numel(Conf)>0
    adjust = true;
else
    adjust = false;
end

load ./MaskenEtc/MNI152GM.mat

mgroup = spm_str_manip(GLMFlags.Group,'rt');

fprintf(1,'%s\n',['Processing RS-data of ' int2str(numel(GLMFlags.sub)) ' subjects ...'  ]);

for ana = 1:numel(GLMFlags.Label)
    if ~(exist(fullfile(cwd,'Projects','RS','GLM',mgroup,GLMFlags.project,des,GLMFlags.Label{ana}))) && ... 
            exist(fullfile(cwd,'Projects','RS','GLM',mgroup,GLMFlags.project,des,[GLMFlags.Label{ana} '.zip']))
        unzip(fullfile(cwd,'Projects','RS','GLM',mgroup,GLMFlags.project,des,[GLMFlags.Label{ana} '.zip']),...
            fullfile(cwd,'Projects','RS','GLM',mgroup,GLMFlags.project,des))
    end
end

for s=1:numel(GLMFlags.sub)
    
    fprintf(1,'%s\n',['Subject ' int2str(s) ' / ' int2str(numel(GLMFlags.sub))])
    
    allDa = 1;
    for ana = 1:numel(GLMFlags.Label)
        if adjust
            if ~(exist(fullfile(cwd,'Projects','RS','GLM',mgroup,GLMFlags.project,des, ... 
                GLMFlags.Label{ana},[GLMFlags.sub{s}  '_raw.nii']))>0)
                allDa = 0;
            end
        else
            if ~(exist(fullfile(cwd,'Projects','RS','GLM',mgroup,GLMFlags.project, ... 
                des,GLMFlags.Label{ana},[GLMFlags.sub{s}  '.nii']))>0)
                allDa = 0;
            end
        end
    end
    
    if allDa == 0
        
        load(fullfile(datadir,'DATA',GLMFlags.SubDir{s},GLMFlags.sub{s},RSdir, ... 
            ['Counfounds_' GLMFlags.sub{s} '.mat']));
        V.fname = fullfile(datadir,'DATA',GLMFlags.SubDir{s},GLMFlags.sub{s},RSdir,spm_str_manip(V.fname,'t'));
        regnow = [reg(:,Qreg) ones(size(reg,1),1)];
        y = spm_read_vols(V);
        
        bb = regnow\y;
        y = y-regnow*bb; clear bb
        y = conn_filter(TR,Filt,y);
        
        for ana = 1:numel(GLMFlags.Label)
            
            [~,~,] = mkdir(fullfile(cwd,'Projects','RS','GLM',mgroup,GLMFlags.project,des,GLMFlags.Label{ana}));
            
            fprintf(1,'%s\n',['...' GLMFlags.Label{ana}])
            
            if ~(exist(fullfile(cwd,'Projects','RS','GLM',mgroup,GLMFlags.project ... 
                ,des,GLMFlags.Label{ana},[GLMFlags.sub{s}  '.nii']))>0)
                
                if GLMFlags.Coordinates==0
                    sXYZ = VM.mat * [maskXYZ; ones(1,size(maskXYZ,2))];
                    Vx = spm_vol(GLMFlags.VOI{ana});
                    D = spm_get_data(Vx,(inv(Vx.mat) * VM.mat * [maskXYZ; ones(1,size(maskXYZ,2))]));
                    ind = find(D>0);
                else
                    sXYZ = VM.mat * [maskXYZ; ones(1,size(maskXYZ,2))];
                    D    =  sqrt((sXYZ(1,:)-GLMFlags.X(ana)).^2+(sXYZ(2,:)-GLMFlags.Y(ana)).^2+(sXYZ(3,:) ... 
                                - GLMFlags.Z(ana)).^2);
                    ind = find(D<GLMFlags.Radius);
                end
                
                switch GLMFlags.gotMask
                    case 1
                        Y3  = ev(y(:,ind));
                    case 2
                        fil = dir(fullfile(datadir,'DATA',GLMFlags.SubDir{s},GLMFlags.sub{s},RSdir,['c1w*.nii']));
                        Vgm = spm_vol(fullfile(datadir,'DATA',GLMFlags.SubDir{s},GLMFlags.sub{s},RSdir,fil(1).name));
                        gm  = spm_get_data(Vgm,round(inv(Vgm(1).mat) * sXYZ(:,ind)));
                        Y3  = ev(y(:,ind(gm>=median(gm))));
                    case 3
                        Vgm = spm_vol(fullfile(cwd,'MaskenEtc','forRS',[spm_str_manip(GLMFlags.Group,'rt') '.nii']));
                        gm  = spm_get_data(Vgm,round(inv(Vgm(1).mat) * sXYZ(:,ind)));
                        Y3  = ev(y(:,ind(gm>=median(gm))));
                end
                
                
                z = nan(1,size(y,2));
                
                Y3  = Y3-mean(Y3);
                Y3  = (Y3./norm(Y3))';
                
                for i=1:size(y,2)
                    z(i) =  Y3 * (y(:,i) ./ norm(y(:,i)));
                end
                
                z = atanh(z);
                z(abs(z)<1E-5) = 1E-5;
                dat = accumarray(maskXYZ',z,VM.dim);
                
                Vo = VM;
                Vo = rmfield(Vo,'pinfo');
                
                if adjust
                    Vo.fname = fullfile(cwd,'Projects','RS','GLM',mgroup,GLMFlags.project, ... 
                        des,GLMFlags.Label{ana},[GLMFlags.sub{s}  '_raw.nii']);
                else
                    Vo.fname = fullfile(cwd,'Projects','RS','GLM',mgroup,GLMFlags.project, ... 
                        des,GLMFlags.Label{ana},[GLMFlags.sub{s}  '.nii']);
                end
                spm_write_vol(Vo,dat);
            end
            
        end
    end
end



if adjust
    OK = 1;
    for s=1:numel(GLMFlags.sub)
        for ana = 1:numel(GLMFlags.Label)
            if ~exist(fullfile(cwd,'Projects','RS','GLM',mgroup,GLMFlags.project,des, ... 
                GLMFlags.Label{ana},[GLMFlags.sub{s}  '.nii']));
                OK=0;
            end
        end
    end
    
    
    Vgm = spm_vol(fullfile(cwd,'MaskenEtc','MNI152GM.nii'));
    gm  = spm_read_vols(Vgm);
    ind = find(gm>0); clear XYZ
    [XYZ(:,1) XYZ(:,2) XYZ(:,3)] = ind2sub(Vgm.dim,ind); clear ind gm
    
    
    
    
    
    if OK == 0
        
        for ana = 1:numel(GLMFlags.Label)
            dat = zeros(numel(GLMFlags.sub),size(XYZ,1));
            for s=1:numel(GLMFlags.sub)
                V        = spm_vol(fullfile(cwd,'Projects','RS','GLM',mgroup,GLMFlags.project, ... 
                    des,GLMFlags.Label{ana},[GLMFlags.sub{s}  '_raw.nii']));
                sXYZ     = V.mat \ Vgm.mat * [XYZ'; ones(1,size(XYZ,1))];
                dat(s,:) = spm_sample_vol(V,XYZ(:,1),XYZ(:,2),XYZ(:,3),1);
            end
            clear sXYZ
            
            
            if GLMFlags.gotConfx>1
                
                mod = mat2cell(GLMFlags.Cov(:,Conf),size(GLMFlags.Cov,1),[ones(1,numel(Conf))]);
                ord = GLMFlags.gotConfx-1;
                useP = 0; try; parfor i=1:10; i; end; useP=1; end
                
                XYZ(sum(isnan(dat))>0,:) = [];
                dat(:,sum(isnan(dat))>0) = [];
                
                c   = mean(dat);
                dat = dat-repmat(c,size(dat,1),1);
                if useP
                    parfor i=1:size(dat,2)
                        [p,table,stats] = ... 
                            anovan(dat(:,i),mod,'model',ord,'display','off','sstype',3,'continuous',cont);
                        dat(:,i) = stats.resid+c(i);
                    end
                else
                    for i=1:size(dat,2)
                        [p,table,stats] = ... 
                            anovan(dat(:,i),mod,'model',ord,'display','off','sstype',3,'continuous',cont);
                        dat(:,i) = stats.resid+c(i);
                    end
                end
                dat = dat+repmat(c,size(dat,1),1);
                
            end
            
            
            % Jetzt wieder rausschreiben
            for s=1:numel(GLMFlags.sub)
                Vo = rmfield(Vgm,'pinfo');
                Vo = rmfield(Vo,'private');
                Vo.fname = fullfile(cwd,'Projects','RS','GLM',mgroup,GLMFlags.project, ... 
                    des,GLMFlags.Label{ana},[GLMFlags.sub{s}  '.nii']);
                Vo        = spm_write_vol(Vo,accumarray(XYZ,dat(s,:),Vo.dim));
            end
            
        end
    end
    clear dat
    
end



if GLMFlags.gotBlock>1
    GLMFlags.xGroup = GLMFlags.Cov(:,GLMFlags.gotBlock-1);
    if ~all(unique(GLMFlags.xGroup)==[1:numel(unique(GLMFlags.xGroup))]')
        tmp = unique(GLMFlags.xGroup);
        vals = unique(tmp);
        GLMFlags.xGroup = ones(1,numel(GLMFlags.xGroup));
        for i=1:numel(vals)
            GLMFlags.xGroup(GLMFlags.xGroup==vals(i)) = i;
        end
    end
else
    GLMFlags.xGroup = ones(size(GLMFlags.Cov,1),1);
end


clear sXYZ

if isfield(GLMFlags,'Label') && numel(GLMFlags.Label)>0 && GLMFlags.gotVBM>1
    
    Vgm = spm_vol(fullfile(cwd,'MaskenEtc','MNI152GM.nii'));
    gm  = spm_read_vols(Vgm);
    ind = find(gm>0); clear XYZ
    [XYZ(:,1) XYZ(:,2) XYZ(:,3)] = ind2sub(Vgm.dim,ind); clear ind gm    
    vois = 1:numel(GLMFlags.Label);
    XYZmm = (Vgm.mat) * [XYZ'; ones(1,size(XYZ,1))];
    
    for ana = 1:numel(vois)
        if GLMFlags.Coordinates == 0
            VOI(ana) = spm_vol(GLMFlags.VOI{vois(ana)});
            dat = find(spm_read_vols(VOI(ana))); clear tmp
            [tmp(:,1) tmp(:,2) tmp(:,3)] = ind2sub(VOI(ana).dim,dat);
            sXYZ{ana} = tmp';
        else
            VOI(ana) = Vgm;
            D = sqrt((XYZmm(1,:)-GLMFlags.X(vois(ana))).^2+(XYZmm(2,:) ... 
                    - GLMFlags.Y(vois(ana))).^2+(XYZmm(3,:)-GLMFlags.Z(vois(ana))).^2);
            sXYZ{ana} = XYZ(D<=3,:)';
        end
    end
    
    dat = zeros(numel(GLMFlags.sub),numel(sXYZ)); do = GLMFlags.gotVBM;
    for s=1:numel(GLMFlags.sub)
        if do == 2
            Vo = spm_vol(fullfile(datadir,'DATA',GLMFlags.SubDir{s},GLMFlags.sub{s},'3D', ... 
                ['mwrp1' GLMFlags.sub{s} '.nii']));
            for ana = 1:numel(sXYZ)
                nowXYZ     = inv(Vo.mat) * VOI(ana).mat * [sXYZ{ana}; ones(1,size(sXYZ{ana},2))];
                dat(s,ana) = sum(spm_get_data(Vo,nowXYZ));
            end
        elseif do == 3
            Vo = spm_vol(fullfile(datadir,'DATA',GLMFlags.SubDir{s},GLMFlags.sub{s},'3D' ... 
                ,['m0wrp1' GLMFlags.sub{s} '.nii']));
            for ana = 1:numel(sXYZ)
                nowXYZ     = inv(Vo.mat) * VOI(ana).mat * [sXYZ{ana}; ones(1,size(sXYZ{ana},2))];
                dat(s,ana) = sum(spm_get_data(Vo,nowXYZ));
            end
        end
    end
    
    for ana = 1:numel(vois)
        GLMFlags.Covariates{size(GLMFlags.Cov,2)+ana} = GLMFlags.Label{ana};
    end
    if isfield(GLMFlags,'ActiveCov')
        GLMFlags.ActiveCov = [GLMFlags.ActiveCov (1:numel(GLMFlags.Label))+size(GLMFlags.Cov,2)];
        GLMFlags.Cov       = [GLMFlags.Cov dat];
    else
        GLMFlags.ActiveCov = [(1:numel(GLMFlags.Label))+size(GLMFlags.Cov,2)];
        GLMFlags.Cov       = [GLMFlags.Cov dat];
        if GLMFlags.gotAnalysis == 1
            GLMFlags.gotAnalysis = 3;
        end
    end
    
end







if any(GLMFlags.isPat)
    analyses = {'GroupDifferences','GroupDifferences_Covariate','Group_x_Covariate',...
        'Patients_Covariate','Patients_Covariates','Paired_t-test','Paired_t-test_Covariate'};
    
    
    if GLMFlags.gotBlock==1
        folderNow = ['ANOVA_' analyses{GLMFlags.gotAnalysis}];
    else
        folderNow = ['ANOVA_' analyses{GLMFlags.gotAnalysis} '_by-' GLMFlags.Covariates{GLMFlags.gotBlock-1}];
    end
    
    
    
    switch GLMFlags.gotAnalysis
        case 1
            load(fullfile(cwd,'matfiles','se_ANOVA.mat'));
            matlabbatch{1}.spm.stats.factorial_design.masking.em{1} = fullfile(cwd,'MaskenEtc','MNI152GM.nii');
            matlabbatch{1}.spm.stats.factorial_design.masking.im = 0;
            
            matlabbatch{1}.spm.stats.factorial_design.dir{1} = fullfile(cwd,'Projects','RS','GLM', ... 
                mgroup,GLMFlags.project,des,folderNow);
            mkdir(matlabbatch{1}.spm.stats.factorial_design.dir{1});
            for s=1:numel(GLMFlags.sub)
                matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(s).scans = {};
                for ana = 1:numel(GLMFlags.Label)
                    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(s).scans{ana,1} = ...
                        fullfile(cwd,'Projects','RS','GLM',mgroup,GLMFlags.project,des,GLMFlags.Label{ana}, ... 
                            [GLMFlags.sub{s}  '.nii']);
                end
                matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(s).conds = ....
                    [1:numel(GLMFlags.Label)]+GLMFlags.isPat(s)*numel(GLMFlags.Label)+(GLMFlags.xGroup(s)-1) ... 
                        * numel(GLMFlags.Label)*2;
            end
             cd(matlabbatch{1}.spm.stats.factorial_design.dir{1}) ; save Design matlabbatch; cd(cwd);
             spm_jobman('run',matlabbatch)
            
            
            load(fullfile(cwd,'matfiles','se_Contrasts.mat'))
            matlabbatch{1}.spm.stats.con.consess = matlabbatch{1}.spm.stats.con.consess(1:2);
            matlabbatch{1}.spm.stats.con.spmmat{1} = fullfile(cwd,'Projects','RS','GLM',mgroup,GLMFlags.project, ... 
                des,folderNow,'SPM.mat');
            matlabbatch{1}.spm.stats.con.consess{1}.fcon.name = 'EoI';
            try; matlabbatch{1}.spm.stats.con.consess{1} =  ... 
                rmfield(matlabbatch{1}.spm.stats.con.consess{1},'tcon'); end
            matlabbatch{1}.spm.stats.con.consess{1}.fcon.convec =  ... 
                {eye(numel(GLMFlags.Label)*2*numel(unique(GLMFlags.xGroup)))};
            matlabbatch{1}.spm.stats.con.consess{2}.fcon.name = 'EoPathology';
            try; matlabbatch{1}.spm.stats.con.consess{2} =  ... 
                rmfield(matlabbatch{1}.spm.stats.con.consess{2},'tcon'); end
            c = zeros(numel(GLMFlags.Label)*numel(unique(GLMFlags.xGroup)), ... 
                numel(GLMFlags.Label)*2*numel(unique(GLMFlags.xGroup)) );
            cnt=1;
            for i= 1:numel(unique(GLMFlags.xGroup))
                for ii= 1:numel(GLMFlags.Label)
                    c(cnt,ii+(i-1)*numel(GLMFlags.Label)*2)=1;    % Controls
                    c(cnt,ii+numel(GLMFlags.Label)+(i-1)*numel(GLMFlags.Label)*2)=-1; % Patients
                    cnt = cnt+1;
                end
            end
            matlabbatch{1}.spm.stats.con.consess{2}.fcon.convec = {c};
            spm_jobman('run',matlabbatch)
            
            load(fullfile(cwd,'matfiles','se_Contrasts.mat')); cnt = 1;
            matlabbatch{1}.spm.stats.con.spmmat{1} =  ... 
                fullfile(cwd,'Projects','RS','GLM',mgroup,GLMFlags.project,des,folderNow,'SPM.mat');
            for ii= 1:numel(GLMFlags.Label)
                matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.name = [GLMFlags.Label{ii} ' +'];
                vec = zeros(1,numel(GLMFlags.Label)*2*numel(unique(GLMFlags.xGroup)));
                for i= 1:numel(unique(GLMFlags.xGroup))
                    vec(ii+(i-1)*numel(GLMFlags.Label)*2)=1;    % Controls
                    vec(ii+numel(GLMFlags.Label)+(i-1)*numel(GLMFlags.Label)*2)=1;% Patients
                end
                matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.convec = vec; cnt = cnt+1;
            end
            for ii= 1:numel(GLMFlags.Label)
                matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.name = [GLMFlags.Label{ii} ' -'];
                vec = zeros(1,numel(GLMFlags.Label)*2*numel(unique(GLMFlags.xGroup)));
                for i= 1:numel(unique(GLMFlags.xGroup))
                    vec(ii+(i-1)*numel(GLMFlags.Label)*2)=-1;    % Controls
                    vec(ii+numel(GLMFlags.Label)+(i-1)*numel(GLMFlags.Label)*2)=-1; % Patients
                end
                matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.convec = vec; cnt = cnt+1;
            end
            for ii= 1:numel(GLMFlags.Label)
                matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.name = [GLMFlags.Label{ii} ' Patients  +'];
                vec = zeros(1,numel(GLMFlags.Label)*2*numel(unique(GLMFlags.xGroup)));
                for i= 1:numel(unique(GLMFlags.xGroup))
                    vec(ii+(i-1)*numel(GLMFlags.Label)*2)=-1;    % Controls
                    vec(ii+numel(GLMFlags.Label)+(i-1)*numel(GLMFlags.Label)*2)=1; % Patients
                end
                matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.convec = vec; cnt = cnt+1;
            end
            for ii= 1:numel(GLMFlags.Label)
                matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.name = [GLMFlags.Label{ii} ' Patients  -'];
                vec = zeros(1,numel(GLMFlags.Label)*2*numel(unique(GLMFlags.xGroup)));
                for i= 1:numel(unique(GLMFlags.xGroup))
                    vec(ii+(i-1)*numel(GLMFlags.Label)*2)=1;    % Controls
                    vec(ii+numel(GLMFlags.Label)+(i-1)*numel(GLMFlags.Label)*2)=-1; % Patients
                end
                matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.convec = vec; cnt = cnt+1;
            end
            spm_jobman('run',matlabbatch)
            
            
            
        case 2
            
            if GLMFlags.gotBlock==1
                folderNow = ['ANOVA_' analyses{GLMFlags.gotAnalysis}];
            else
                folderNow = ['ANOVA_' analyses{GLMFlags.gotAnalysis} '_by-' ... 
                    GLMFlags.Covariates{GLMFlags.gotBlock-1}];
            end
            
            load(fullfile(cwd,'matfiles','se_ANOVAcov.mat'))
            matlabbatch{1}.spm.stats.factorial_design.masking.em{1} = fullfile(cwd,'MaskenEtc','MNI152GM.nii');
            matlabbatch{1}.spm.stats.factorial_design.masking.im = 0;
            matlabbatch{1}.spm.stats.factorial_design.dir{1} =  ... 
                fullfile(cwd,'Projects','RS','GLM',mgroup,GLMFlags.project,des,folderNow);
            mkdir(matlabbatch{1}.spm.stats.factorial_design.dir{1});
            for s=1:numel(GLMFlags.sub)
                matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(s).scans = {};
                for ana = 1:numel(GLMFlags.Label)
                    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(s).scans{ana,1} = ...
                        fullfile(cwd,'Projects','RS','GLM',mgroup,GLMFlags.project,des,GLMFlags.Label{ana}, ... 
                        [GLMFlags.sub{s} '.nii']);
                end
                matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(s).conds = ....
                    [1:numel(GLMFlags.Label)]+GLMFlags.isPat(s)*numel(GLMFlags.Label)+(GLMFlags.xGroup(s)-1) ... 
                    * numel(GLMFlags.Label)*2;
            end
            
            for cx = 1:numel(GLMFlags.ActiveCov)
                matlabbatch{1}.spm.stats.factorial_design.cov(cx) = ... 
                    matlabbatch{1}.spm.stats.factorial_design.cov(1);
                c = GLMFlags.Cov(:,GLMFlags.ActiveCov(cx));
                c(isnan(c) & GLMFlags.isPat'==0) = mean(c(~isnan(c) & GLMFlags.isPat'==0));
                c(isnan(c) & GLMFlags.isPat'==1) = mean(c(~isnan(c) & GLMFlags.isPat'==1));
                c(GLMFlags.isPat'==0) = c(GLMFlags.isPat'==0)-mean(c(GLMFlags.isPat'==0));
                c(GLMFlags.isPat'==1) = c(GLMFlags.isPat'==1)-mean(c(GLMFlags.isPat'==1));
                c = repmat(c',numel(GLMFlags.Label),1);
                matlabbatch{1}.spm.stats.factorial_design.cov(cx).c = c (:);
                matlabbatch{1}.spm.stats.factorial_design.cov(cx).cname =  ... 
                    GLMFlags.Covariates{GLMFlags.ActiveCov(cx)};
            end
             cd(matlabbatch{1}.spm.stats.factorial_design.dir{1}) ; save Design matlabbatch; cd(cwd);
            save Design matlabbatch
            spm_jobman('run',matlabbatch)
            
            
            
            load(fullfile(cwd,'matfiles','se_Contrasts.mat'))
            matlabbatch{1}.spm.stats.con.consess = matlabbatch{1}.spm.stats.con.consess(1:2);
            matlabbatch{1}.spm.stats.con.spmmat{1} =  ... 
                fullfile(cwd,'Projects','RS','GLM',mgroup,GLMFlags.project,des,folderNow,'SPM.mat');
            matlabbatch{1}.spm.stats.con.consess{1}.fcon.name = 'EoI';
            try; matlabbatch{1}.spm.stats.con.consess{1} =  ... 
                rmfield(matlabbatch{1}.spm.stats.con.consess{1},'tcon'); end
            matlabbatch{1}.spm.stats.con.consess{1}.fcon.convec =  ... 
                {eye(numel(GLMFlags.Label)*2*numel(unique(GLMFlags.xGroup)))};
            matlabbatch{1}.spm.stats.con.consess{2}.fcon.name = 'EoPathology';
            try; matlabbatch{1}.spm.stats.con.consess{2} =  ... 
                rmfield(matlabbatch{1}.spm.stats.con.consess{2},'tcon'); end
            c = zeros(numel(GLMFlags.Label)*numel(unique(GLMFlags.xGroup)), ... 
                numel(GLMFlags.Label)*2*numel(unique(GLMFlags.xGroup)) );
            cnt=1;
            for i= 1:numel(unique(GLMFlags.xGroup))
                for ii= 1:numel(GLMFlags.Label)
                    c(cnt,ii+(i-1)*numel(GLMFlags.Label)*2)=1;    % Controls
                    c(cnt,ii+numel(GLMFlags.Label)+(i-1)*numel(GLMFlags.Label)*2)=-1; % Patients
                    cnt = cnt+1;
                end
            end
            matlabbatch{1}.spm.stats.con.consess{2}.fcon.convec = {c}; %%%% matlabbatch{1}.spm.stats.con.delete = 1;
            spm_jobman('run',matlabbatch)
            
            
            load(fullfile(cwd,'matfiles','se_Contrasts.mat')); cnt = 1;
            matlabbatch{1}.spm.stats.con.spmmat{1} =  ... 
                fullfile(cwd,'Projects','RS','GLM',mgroup,GLMFlags.project,des,folderNow,'SPM.mat');
            for ii= 1:numel(GLMFlags.Label)
                matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.name = [GLMFlags.Label{ii} ' +'];
                vec = zeros(1,numel(GLMFlags.Label)*2*numel(unique(GLMFlags.xGroup)));
                for i= 1:numel(unique(GLMFlags.xGroup))
                    vec(ii+(i-1)*numel(GLMFlags.Label)*2)=1;    % Controls
                    vec(ii+numel(GLMFlags.Label)+(i-1)*numel(GLMFlags.Label)*2)=1;% Patients
                end
                matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.convec = vec; cnt = cnt+1;
            end
            for ii= 1:numel(GLMFlags.Label)
                matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.name = [GLMFlags.Label{ii} ' -'];
                vec = zeros(1,numel(GLMFlags.Label)*2*numel(unique(GLMFlags.xGroup)));
                for i= 1:numel(unique(GLMFlags.xGroup))
                    vec(ii+(i-1)*numel(GLMFlags.Label)*2)=-1;    % Controls
                    vec(ii+numel(GLMFlags.Label)+(i-1)*numel(GLMFlags.Label)*2)=-1; % Patients
                end
                matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.convec = vec; cnt = cnt+1;
            end
            for ii= 1:numel(GLMFlags.Label)
                matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.name = [GLMFlags.Label{ii} ' Patients  +'];
                vec = zeros(1,numel(GLMFlags.Label)*2*numel(unique(GLMFlags.xGroup)));
                for i= 1:numel(unique(GLMFlags.xGroup))
                    vec(ii+(i-1)*numel(GLMFlags.Label)*2)=-1;    % Controls
                    vec(ii+numel(GLMFlags.Label)+(i-1)*numel(GLMFlags.Label)*2)=1; % Patients
                end
                matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.convec = vec; cnt = cnt+1;
            end
            for ii= 1:numel(GLMFlags.Label)
                matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.name = [GLMFlags.Label{ii} ' Patients  -'];
                vec = zeros(1,numel(GLMFlags.Label)*2*numel(unique(GLMFlags.xGroup)));
                for i= 1:numel(unique(GLMFlags.xGroup))
                    vec(ii+(i-1)*numel(GLMFlags.Label)*2)=1;    % Controls
                    vec(ii+numel(GLMFlags.Label)+(i-1)*numel(GLMFlags.Label)*2)=-1; % Patients
                end
                matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.convec = vec; cnt = cnt+1;
            end
            spm_jobman('run',matlabbatch)
            
            load(fullfile(cwd,'matfiles','se_Contrasts.mat'))
            matlabbatch{1}.spm.stats.con.spmmat{1} =  ... 
                fullfile(cwd,'Projects','RS','GLM',mgroup,GLMFlags.project,des,folderNow,'SPM.mat'); cnt = 1;
            for cx = 1:numel(GLMFlags.ActiveCov)
                for ana = 1:numel(GLMFlags.Label)
                    matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.name =  ... 
                        [GLMFlags.Label{ana} ' ' GLMFlags.Covariates{GLMFlags.ActiveCov(cx)} ' +'];
                    vec = zeros(1,numel(GLMFlags.Label)*(2+numel(GLMFlags.ActiveCov)));
                    vec(ana+(numel(GLMFlags.Label)*(1+cx)))  = 1;
                    matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.convec = vec;
                    cnt = cnt+1;
                    matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.name =  ... 
                        [GLMFlags.Label{ana} ' ' GLMFlags.Covariates{GLMFlags.ActiveCov(cx)} ' -'];
                    vec = zeros(1,numel(GLMFlags.Label)*(2+numel(GLMFlags.ActiveCov)));
                    vec(ana+(numel(GLMFlags.Label)*(1+cx)))  = -1;
                    matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.convec = vec;
                    cnt = cnt+1;
                    matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.name =  ... 
                        [GLMFlags.Label{ana} ' ' GLMFlags.Covariates{GLMFlags.ActiveCov(cx)} ' Pat+'];
                    vec = zeros(1,numel(GLMFlags.Label)*(2+numel(GLMFlags.ActiveCov)));
                    vec(ana+(numel(GLMFlags.Label)*(2+cx)))  = 1;
                    matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.convec = vec;
                    cnt = cnt+1;
                    matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.name =  ... 
                        [GLMFlags.Label{ana} ' ' GLMFlags.Covariates{GLMFlags.ActiveCov(cx)} ' Pat-'];
                    vec = zeros(1,numel(GLMFlags.Label)*(2+numel(GLMFlags.ActiveCov)));
                    vec(ana+(numel(GLMFlags.Label)*(2+cx)))  = -1;
                    matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.convec = vec;
                    cnt = cnt+1;
                end
            end
            spm_jobman('run',matlabbatch)
            
            
            
        case 3
            load(fullfile(cwd,'matfiles','se_ANOVAcov.mat'))
            matlabbatch{1}.spm.stats.factorial_design.masking.em{1} = fullfile(cwd,'MaskenEtc','MNI152GM.nii');
            matlabbatch{1}.spm.stats.factorial_design.masking.im = 0;
            matlabbatch{1}.spm.stats.factorial_design.dir{1} =  ... 
                fullfile(cwd,'Projects','RS','GLM',mgroup,GLMFlags.project,des,folderNow);
            mkdir(matlabbatch{1}.spm.stats.factorial_design.dir{1});
            for s=1:numel(GLMFlags.sub)
                matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(s).scans = {};
                for ana = 1:numel(GLMFlags.Label)
                    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(s).scans{ana,1} = ...
                        fullfile(cwd,'Projects','RS','GLM',mgroup,GLMFlags.project, ... 
                        des,GLMFlags.Label{ana},[GLMFlags.sub{s}  '.nii']);
                end
                matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(s).conds =  ... 
                [1:numel(GLMFlags.Label)]+(GLMFlags.isPat(s)*numel(GLMFlags.Label));
            end
            for cx = 1:numel(GLMFlags.ActiveCov)
                matlabbatch{1}.spm.stats.factorial_design.cov(cx) =  ... 
                    matlabbatch{1}.spm.stats.factorial_design.cov(1);
                c = GLMFlags.Cov(:,GLMFlags.ActiveCov(cx));
                c(isnan(c) & GLMFlags.isPat'==0) = mean(c(~isnan(c) & GLMFlags.isPat'==0));
                c(isnan(c) & GLMFlags.isPat'==1) = mean(c(~isnan(c) & GLMFlags.isPat'==1));
                c = repmat(c',numel(GLMFlags.Label),1);
                matlabbatch{1}.spm.stats.factorial_design.cov(cx).c = c (:);
                matlabbatch{1}.spm.stats.factorial_design.cov(cx).cname =  ... 
                GLMFlags.Covariates{GLMFlags.ActiveCov(cx)};
            end
            cd(matlabbatch{1}.spm.stats.factorial_design.dir{1});
            save Design matlabbatch; cd(cwd);
            spm_jobman('run',matlabbatch)
            
            
            load(fullfile(cwd,'matfiles','se_Contrasts.mat')); cnt = 1;
            matlabbatch{1}.spm.stats.con.spmmat{1} =  ... 
                fullfile(cwd,'Projects','RS','GLM',mgroup,GLMFlags.project,des,folderNow,'SPM.mat');
            for ana = 1:numel(GLMFlags.Label)
                matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.name = [GLMFlags.Label{ana} ' Con +'];
                vec = zeros(1,numel(GLMFlags.Label)*2+numel(GLMFlags.Label)*2*numel(GLMFlags.ActiveCov));
                vec(ana)  = 1;  matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.convec = vec; cnt = cnt+1;
                matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.name = [GLMFlags.Label{ana} ' Con -'];
                vec = zeros(1,numel(GLMFlags.Label)*2+numel(GLMFlags.Label)*2*numel(GLMFlags.ActiveCov));
                vec(ana)  = -1; matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.convec = vec; cnt = cnt+1;
                matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.name = [GLMFlags.Label{ana} ' Pat +'];
                vec = zeros(1,numel(GLMFlags.Label)*2+numel(GLMFlags.Label)*2*numel(GLMFlags.ActiveCov));
                vec(ana+numel(GLMFlags.Label))  = 1; 
                matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.convec = vec;
                cnt = cnt+1;
                matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.name = [GLMFlags.Label{ana} ' Pat -'];
                vec = zeros(1,numel(GLMFlags.Label)*2+numel(GLMFlags.Label)*2*numel(GLMFlags.ActiveCov));
                vec(ana+numel(GLMFlags.Label))  = -1;
                matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.convec = vec;
                cnt = cnt+1;
            end
            spm_jobman('run',matlabbatch)
            
            load(fullfile(cwd,'matfiles','se_Contrasts.mat')); cnt = 1;
            matlabbatch{1}.spm.stats.con.spmmat{1} =  ... 
                fullfile(cwd,'Projects','RS','GLM',mgroup,GLMFlags.project,des,folderNow,'SPM.mat');
            for ana = 1:numel(GLMFlags.Label)
                matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.name = [GLMFlags.Label{ana} ' Pat > Con'];
                vec = zeros(1,numel(GLMFlags.Label)*2+numel(GLMFlags.Label)*2*numel(GLMFlags.ActiveCov));
                vec(ana)  = -1; vec(ana+numel(GLMFlags.Label))  = 1;
                matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.convec = vec;
                cnt = cnt+1;
                matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.name = [GLMFlags.Label{ana} ' Pat < Con'];
                vec = zeros(1,numel(GLMFlags.Label)*2+numel(GLMFlags.Label)*2*numel(GLMFlags.ActiveCov));
                vec(ana)  = 1; vec(ana+numel(GLMFlags.Label))  = -1;
                matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.convec = vec;
                cnt = cnt+1;
            end
            spm_jobman('run',matlabbatch)
            
            load(fullfile(cwd,'matfiles','se_Contrasts.mat')); cnt = 1;
            matlabbatch{1}.spm.stats.con.spmmat{1} =  ... 
                fullfile(cwd,'Projects','RS','GLM',mgroup,GLMFlags.project,des,folderNow,'SPM.mat');
            for cx = 1:numel(GLMFlags.ActiveCov)
                for ana = 1:numel(GLMFlags.Label)
                    matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.name =  ... 
                        [GLMFlags.Label{ana} ' ' GLMFlags.Covariates{GLMFlags.ActiveCov(cx)} ' Con +'];
                    vec = zeros(1,numel(GLMFlags.Label)*2+numel(GLMFlags.Label)*2*numel(GLMFlags.ActiveCov));
                    vec(ana+(numel(GLMFlags.Label)*(1+cx)))  = 1;
                    matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.convec = vec;
                    cnt = cnt+1;
                    matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.name =  ... 
                        [GLMFlags.Label{ana} ' ' GLMFlags.Covariates{GLMFlags.ActiveCov(cx)} ' Con -'];
                    vec = zeros(1,numel(GLMFlags.Label)*2+numel(GLMFlags.Label)*2*numel(GLMFlags.ActiveCov));
                    vec(ana+(numel(GLMFlags.Label)*(1+cx))) = -1;
                    matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.convec = vec;
                    cnt = cnt+1;
                    matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.name =  ... 
                        [GLMFlags.Label{ana} ' ' GLMFlags.Covariates{GLMFlags.ActiveCov(cx)} ' Pat +'];
                    vec = zeros(1,numel(GLMFlags.Label)*2+numel(GLMFlags.Label)*2*numel(GLMFlags.ActiveCov));
                    vec(ana+numel(GLMFlags.Label)+(numel(GLMFlags.Label)*(1+cx)))  = 1;
                    matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.convec = vec;
                    cnt = cnt+1;
                    matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.name = ... 
                        [GLMFlags.Label{ana} ' ' GLMFlags.Covariates{GLMFlags.ActiveCov(cx)} ' Pat -'];
                    vec = zeros(1,numel(GLMFlags.Label)*2+numel(GLMFlags.Label)*2*numel(GLMFlags.ActiveCov));
                    vec(ana+numel(GLMFlags.Label)+(numel(GLMFlags.Label)*(1+cx))) = -1;
                    matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.convec = vec;
                    cnt = cnt+1;
                end
                spm_jobman('run',matlabbatch)
            end
            
            load(fullfile(cwd,'matfiles','se_Contrasts.mat')); cnt = 1;
            matlabbatch{1}.spm.stats.con.spmmat{1} =  ... 
                fullfile(cwd,'Projects','RS','GLM',mgroup,GLMFlags.project,des,folderNow,'SPM.mat');
            for ana = 1:numel(GLMFlags.Label)
                matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.name =  ... 
                [GLMFlags.Label{ana} ' ' GLMFlags.Covariates{GLMFlags.ActiveCov(cx)} ' Pat > Con'];
                vec = zeros(1,numel(GLMFlags.Label)*2+numel(GLMFlags.Label)*2*numel(GLMFlags.ActiveCov));
                vec(ana+(numel(GLMFlags.Label)*(1+cx)))  = -1;
                vec(ana+numel(GLMFlags.Label)+(numel(GLMFlags.Label)*(1+cx)))  = 1;
                matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.convec = vec;
                cnt = cnt+1;
                matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.name = ... 
                    [GLMFlags.Label{ana} ' ' GLMFlags.Covariates{GLMFlags.ActiveCov(cx)} ' Pat < Con'];
                vec = zeros(1,numel(GLMFlags.Label)*2+numel(GLMFlags.Label)*2*numel(GLMFlags.ActiveCov));
                vec(ana+(numel(GLMFlags.Label)*(1+cx)))  = 1;
                vec(ana+numel(GLMFlags.Label)+(numel(GLMFlags.Label)*(1+cx)))  = -1;
                matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.convec = vec;
                cnt = cnt+1;
            end
            spm_jobman('run',matlabbatch)
            
            
            
        case 4
            xfolder = folderNow;
            
            for cx = 1:numel(GLMFlags.ActiveCov)
                load(fullfile(cwd,'matfiles','se_ANOVAcov.mat'))
                matlabbatch{1}.spm.stats.factorial_design.masking.em{1} = fullfile(cwd,'MaskenEtc','MNI152GM.nii');
                matlabbatch{1}.spm.stats.factorial_design.masking.im = 0;
                folderNow = [xfolder '_' GLMFlags.Covariates{GLMFlags.ActiveCov(cx)}];
                matlabbatch{1}.spm.stats.factorial_design.dir{1} = ... 
                    fullfile(cwd,'Projects','RS','GLM',mgroup,GLMFlags.project,des,folderNow);
                mkdir(matlabbatch{1}.spm.stats.factorial_design.dir{1});
                for s=1:numel(GLMFlags.sub)
                    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(s).scans = {};
                    for ana = 1:numel(GLMFlags.Label)
                        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(s).scans{ana,1} = ...
                            fullfile(cwd,'Projects','RS','GLM',mgroup,GLMFlags.project, ... 
                            des,GLMFlags.Label{ana},[GLMFlags.sub{s}  '.nii']);
                    end
                    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(s).conds = ... 
                        [1:numel(GLMFlags.Label)];
                end
                matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject = ...
                    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(GLMFlags.isPat>0);
                matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject = ...
                    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject ... 
                    (find(sum(isnan(GLMFlags.Cov(GLMFlags.isPat>0,GLMFlags.ActiveCov(cx))),2)==0));
                
                c = GLMFlags.Cov(:,GLMFlags.ActiveCov(cx));
                c =c(GLMFlags.isPat>0);
                c =c(find(sum(isnan(GLMFlags.Cov(GLMFlags.isPat>0,GLMFlags.ActiveCov(cx))),2)==0));
                c = repmat(c',numel(GLMFlags.Label),1);
                matlabbatch{1}.spm.stats.factorial_design.cov(1).c = c (:);
                matlabbatch{1}.spm.stats.factorial_design.cov(1).cname = ... 
                GLMFlags.Covariates{GLMFlags.ActiveCov(cx)};
                cd(matlabbatch{1}.spm.stats.factorial_design.dir{1});
                save Design matlabbatch; cd(cwd);
                spm_jobman('run',matlabbatch)
                
                
                load(fullfile(cwd,'matfiles','se_Contrasts.mat'))
                matlabbatch{1}.spm.stats.con.spmmat{1} =  ... 
                    fullfile(cwd,'Projects','RS','GLM',mgroup,GLMFlags.project,des,folderNow,'SPM.mat');
                matlabbatch{1}.spm.stats.con.consess{1}.fcon.name = 'EoI';
                try; matlabbatch{1}.spm.stats.con.consess{1} =  ... 
                    rmfield(matlabbatch{1}.spm.stats.con.consess{1},'tcon'); end
                matlabbatch{1}.spm.stats.con.consess{1}.fcon.convec = {eye(numel(GLMFlags.Label)*2)};
                matlabbatch{1}.spm.stats.con.consess{2}.fcon.name = 'EoConnectivity';
                try; matlabbatch{1}.spm.stats.con.consess{2} = ... 
                    rmfield(matlabbatch{1}.spm.stats.con.consess{2},'tcon'); end
                matlabbatch{1}.spm.stats.con.consess{2}.fcon.convec =  ... 
                    {[eye(numel(GLMFlags.Label)) zeros(numel(GLMFlags.Label))]};
                matlabbatch{1}.spm.stats.con.consess{3}.fcon.name = 'EoCovariate';
                try; matlabbatch{1}.spm.stats.con.consess{2} =  ... 
                    rmfield(matlabbatch{1}.spm.stats.con.consess{3},'tcon'); end
                matlabbatch{1}.spm.stats.con.consess{3}.fcon.convec =  ... 
                    {[zeros(numel(GLMFlags.Label)) eye(numel(GLMFlags.Label))]};
                spm_jobman('run',matlabbatch)
                
                load(fullfile(cwd,'matfiles','se_Contrasts.mat'))
                matlabbatch{1}.spm.stats.con.spmmat{1} =  ... 
                    fullfile(cwd,'Projects','RS','GLM',mgroup,GLMFlags.project,des,folderNow,'SPM.mat');
                for ana = 1:numel(GLMFlags.Label)
                    matlabbatch{1}.spm.stats.con.consess{ana}.tcon.name = [GLMFlags.Label{ana} ' +'];
                    vec = zeros(1,numel(GLMFlags.Label)*2); vec(ana)  = 1;
                    matlabbatch{1}.spm.stats.con.consess{ana}.tcon.convec = vec;
                end
                for ana = 1:numel(GLMFlags.Label)
                    matlabbatch{1}.spm.stats.con.consess{ana+(1*numel(GLMFlags.Label))}.tcon.name =  ... 
                        [GLMFlags.Label{ana} ' -'];
                    vec = zeros(1,numel(GLMFlags.Label)*2); vec(ana)  = -1;
                    matlabbatch{1}.spm.stats.con.consess{ana+(1*numel(GLMFlags.Label))}.tcon.convec = vec;
                end
                for ana = 1:numel(GLMFlags.Label)
                    matlabbatch{1}.spm.stats.con.consess{ana+(2*numel(GLMFlags.Label))}.tcon.name =  ... 
                        [GLMFlags.Label{ana} ': ' GLMFlags.Covariates{GLMFlags.ActiveCov(cx)} ' +'];
                    vec = zeros(1,numel(GLMFlags.Label)*2); vec(ana+numel(GLMFlags.Label))  = 1;
                    matlabbatch{1}.spm.stats.con.consess{ana+(2*numel(GLMFlags.Label))}.tcon.convec = vec;
                end
                for ana = 1:numel(GLMFlags.Label)
                    matlabbatch{1}.spm.stats.con.consess{ana+(3*numel(GLMFlags.Label))}.tcon.name =  ... 
                        [GLMFlags.Label{ana} ': ' GLMFlags.Covariates{GLMFlags.ActiveCov(cx)} ' -'];
                    vec = zeros(1,numel(GLMFlags.Label)*2); vec(ana+numel(GLMFlags.Label))  = -1;
                    matlabbatch{1}.spm.stats.con.consess{ana+(3*numel(GLMFlags.Label))}.tcon.convec = vec;
                end
                spm_jobman('run',matlabbatch)
                
            end
            
            
            
        case 5
            load(fullfile(cwd,'matfiles','se_ANOVAcov.mat'))
            matlabbatch{1}.spm.stats.factorial_design.masking.em{1} = fullfile(cwd,'MaskenEtc','MNI152GM.nii');
            matlabbatch{1}.spm.stats.factorial_design.masking.im = 0;
            matlabbatch{1}.spm.stats.factorial_design.dir{1} =  ... 
                fullfile(cwd,'Projects','RS','GLM',mgroup,GLMFlags.project,des,folderNow);
            mkdir(matlabbatch{1}.spm.stats.factorial_design.dir{1});
            for s=1:numel(GLMFlags.sub)
                matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(s).scans = {};
                for ana = 1:numel(GLMFlags.Label)
                    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(s).scans{ana,1} = ...
                        fullfile(cwd,'Projects','RS','GLM',mgroup,GLMFlags.project, ... 
                        des,GLMFlags.Label{ana},[GLMFlags.sub{s}  '.nii']);
                end
                matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(s).conds =  ... 
                    [1:numel(GLMFlags.Label)];
            end
            matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject =...
                matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(GLMFlags.isPat>0);
            matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject =...
                matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject ... 
                (find(sum(isnan(GLMFlags.Cov(GLMFlags.isPat>0,GLMFlags.ActiveCov)),2)==0));
            
            for cx = 1:numel(GLMFlags.ActiveCov)
                matlabbatch{1}.spm.stats.factorial_design.cov(cx) = ... 
                    matlabbatch{1}.spm.stats.factorial_design.cov(1);
                c = GLMFlags.Cov(:,GLMFlags.ActiveCov(cx));
                c =c(GLMFlags.isPat>0);
                c =c(find(sum(isnan(GLMFlags.Cov(GLMFlags.isPat>0,GLMFlags.ActiveCov)),2)==0));
                c = repmat(c',numel(GLMFlags.Label),1);
                matlabbatch{1}.spm.stats.factorial_design.cov(cx).c = c (:);
                matlabbatch{1}.spm.stats.factorial_design.cov(cx).cname =  ... 
                    GLMFlags.Covariates{GLMFlags.ActiveCov(cx)};
            end
            %      matlabbatch = se_adjustAnova(matlabbatch);
             cd(matlabbatch{1}.spm.stats.factorial_design.dir{1}) ; save Design matlabbatch; cd(cwd);
             spm_jobman('run',matlabbatch)
            
            
            load(fullfile(cwd,'matfiles','se_Contrasts.mat')); cnt = 1;
            matlabbatch{1}.spm.stats.con.spmmat{1} = ... 
                    fullfile(cwd,'Projects','RS','GLM',mgroup,GLMFlags.project,des,folderNow,'SPM.mat');
            for ana = 1:numel(GLMFlags.Label)
                matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.name = [GLMFlags.Label{ana} ' +'];
                vec = zeros(1,numel(GLMFlags.Label)*(numel(GLMFlags.ActiveCov)+1));
                vec(ana)  = 1;  matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.convec = vec; cnt = cnt+1;
                matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.name = [GLMFlags.Label{ana} ' -'];
                vec = zeros(1,numel(GLMFlags.Label)*(numel(GLMFlags.ActiveCov)+1));
                vec(ana)  = -1; matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.convec = vec; cnt = cnt+1;
            end
            spm_jobman('run',matlabbatch)
            
            load(fullfile(cwd,'matfiles','se_Contrasts.mat'))
            matlabbatch{1}.spm.stats.con.spmmat{1} =   ... 
                fullfile(cwd,'Projects','RS','GLM',mgroup,GLMFlags.project,des,folderNow,'SPM.mat');
            cnt = 1;
            for ana = 1:numel(GLMFlags.Label)
                for i=1:numel(GLMFlags.ActiveCov)
                    matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.name =   ... 
                        [GLMFlags.Label{ana} ': ' GLMFlags.Covariates{GLMFlags.ActiveCov(i)} ' +'];
                    vec = zeros(1,numel(GLMFlags.Label)*(numel(GLMFlags.ActiveCov)+1));
                    vec(((i)*numel(GLMFlags.Label))+ana)  = 1;
                    matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.convec = vec;
                    cnt = cnt+1;
                    
                    matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.name =   ... 
                        [GLMFlags.Label{ana} ': ' GLMFlags.Covariates{GLMFlags.ActiveCov(i)}  ' -'];
                    vec = zeros(1,numel(GLMFlags.Label)*(numel(GLMFlags.ActiveCov)+1));
                    vec(((i)*numel(GLMFlags.Label))+ana)  = -1;
                    matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.convec = vec;
                    cnt = cnt+1;
                end
            end
            spm_jobman('run',matlabbatch)
            
            load(fullfile(cwd,'matfiles','se_Contrasts.mat'))
            matlabbatch{1}.spm.stats.con.spmmat{1} =   ... 
                fullfile(cwd,'Projects','RS','GLM',mgroup,GLMFlags.project,des,folderNow,'SPM.mat');
            matlabbatch{1}.spm.stats.con.consess{1}.fcon.name = 'EoConnectivity';
            try; matlabbatch{1}.spm.stats.con.consess{1} =   ... 
                rmfield(matlabbatch{1}.spm.stats.con.consess{1},'tcon'); end
            matlabbatch{1}.spm.stats.con.consess{1}.fcon.convec =   ... 
                {[eye(numel(GLMFlags.Label)) repmat(zeros(numel(GLMFlags.Label)),1,numel(GLMFlags.ActiveCov))]};
            matlabbatch{1}.spm.stats.con.consess{2}.fcon.name = 'EoCovariates';
            try; matlabbatch{1}.spm.stats.con.consess{2} =   ... 
                rmfield(matlabbatch{1}.spm.stats.con.consess{2},'tcon'); end
            matlabbatch{1}.spm.stats.con.consess{2}.fcon.convec =   ... 
                    {[zeros(numel(GLMFlags.Label)*numel(GLMFlags.ActiveCov),numel(GLMFlags.Label)) ... 
                        eye(numel(GLMFlags.Label)*numel(GLMFlags.ActiveCov))]};
            
            cnt = 3;
            for cx=1:numel(GLMFlags.ActiveCov)
                matlabbatch{1}.spm.stats.con.consess{cnt}.fcon.name =   ... 
                    ['Eo' GLMFlags.Covariates{GLMFlags.ActiveCov(cx)}];
                try; matlabbatch{1}.spm.stats.con.consess{cnt} =   ... 
                    rmfield(matlabbatch{1}.spm.stats.con.consess{cnt},'tcon'); end
                vec = zeros(numel(GLMFlags.Label),numel(GLMFlags.Label)*(numel(GLMFlags.ActiveCov)+1));
                for ana = 1:numel(GLMFlags.Label); vec(ana,((cx)*numel(GLMFlags.Label))+ana)  = 1; end
                matlabbatch{1}.spm.stats.con.consess{cnt}.fcon.convec = {vec};
                cnt = cnt+1;
            end
            spm_jobman('run',matlabbatch)
  
        case 6
            for ana = 1:numel(GLMFlags.Label)

                load(fullfile(cwd,'matfiles','Pairedttest.mat'));
                matlabbatch{1}.spm.stats.factorial_design.masking.em{1} = fullfile(cwd,'MaskenEtc','MNI152GM.nii');
                matlabbatch{1}.spm.stats.factorial_design.masking.im = 0;

                matlabbatch{1}.spm.stats.factorial_design.dir{1} = ... 
                    fullfile(cwd,'Projects','RS','GLM',mgroup,GLMFlags.project,des,folderNow,GLMFlags.Label{ana});
                mkdir(matlabbatch{1}.spm.stats.factorial_design.dir{1});
                for s=1:numel(GLMFlags.sub)/2
                    matlabbatch{1}.spm.stats.factorial_design.des.pt.pair(s).scans = {};
                    for p = 1:2
                        matlabbatch{1}.spm.stats.factorial_design.des.pt.pair(s).scans{p,1} = ...
                            fullfile(cwd,'Projects','RS','GLM',mgroup,GLMFlags.project, ... 
                            des,GLMFlags.Label{ana},[GLMFlags.sub{GLMFlags.matchQ(s,p)} '.nii']); 
                    end
                end

                 cd(matlabbatch{1}.spm.stats.factorial_design.dir{1}) ; save Design matlabbatch; cd(cwd);
                spm_jobman('run',matlabbatch)
            end    
                
            for ana = 1:numel(GLMFlags.Label)
                load(fullfile(cwd,'matfiles','se_Contrasts.mat'))
                matlabbatch{1}.spm.stats.con.spmmat{1} =  ... 
                    fullfile(cwd,'Projects','RS','GLM',mgroup,GLMFlags.project, ... 
                    des,folderNow,GLMFlags.Label{ana},'SPM.mat');
                matlabbatch{1}.spm.stats.con.delete = 1;
                
                matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'T1-T2';
                matlabbatch{1}.spm.stats.con.consess{1}.tcon.convec = [1 -1];

                matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'T2-T1';
                matlabbatch{1}.spm.stats.con.consess{2}.tcon.convec = [-1 1];

                prs=numel(GLMFlags.sub)/2;
                matlabbatch{1}.spm.stats.con.consess{3}.tcon.name = 'main effect T1+';
                matlabbatch{1}.spm.stats.con.consess{3}.tcon.convec = [1 0 ones(1,prs)/prs];
              
                matlabbatch{1}.spm.stats.con.consess{4}.tcon.name = 'main effect T1-';
                matlabbatch{1}.spm.stats.con.consess{4}.tcon.convec = [-1 0 -ones(1,prs)/prs];

                matlabbatch{1}.spm.stats.con.consess{5}.tcon.name = 'main effect T2+';
                matlabbatch{1}.spm.stats.con.consess{5}.tcon.convec = [1 0 ones(1,prs)/prs];

                matlabbatch{1}.spm.stats.con.consess{6}.tcon.name = 'main effect T2-';
                matlabbatch{1}.spm.stats.con.consess{6}.tcon.convec = [-1 0 -ones(1,prs)/prs];

                spm_jobman('run',matlabbatch)
            end

            
            
        case 7
            for ana = 1:numel(GLMFlags.Label)

                load(fullfile(cwd,'matfiles','Pairedttestcov.mat'));
                matlabbatch{1}.spm.stats.factorial_design.masking.em{1} = fullfile(cwd,'MaskenEtc','MNI152GM.nii');
                matlabbatch{1}.spm.stats.factorial_design.masking.im = 0;

                matlabbatch{1}.spm.stats.factorial_design.dir{1} = ... 
                    fullfile(cwd,'Projects','RS','GLM',mgroup,GLMFlags.project,des,folderNow,[GLMFlags.Label{ana}]);
                mkdir(matlabbatch{1}.spm.stats.factorial_design.dir{1});
                for s=1:numel(GLMFlags.sub)/2
                    matlabbatch{1}.spm.stats.factorial_design.des.pt.pair(s).scans = {};
                    for p = 1:2
                        matlabbatch{1}.spm.stats.factorial_design.des.pt.pair(s).scans{p,1} = ...
                            fullfile(cwd,'Projects','RS','GLM',mgroup,GLMFlags.project,... 
                            des,GLMFlags.Label{ana},[GLMFlags.sub{GLMFlags.matchQ(s,p)} '.nii']); 
                    end
                end
                
                succ=GLMFlags.matchQ';
                
                for cx = 1:numel(GLMFlags.ActiveCov)
                    c(1,:) = GLMFlags.Cov(succ(1,:),GLMFlags.ActiveCov(cx));
                    c(2,:) = GLMFlags.Cov(succ(2,:),GLMFlags.ActiveCov(cx));
                    if sum(isnan(c(1,:)))==numel(c(1,:))
                         c(1,find(isnan(c(1,:))))=0;
                    else
                         c(1,find(isnan(c(1,:))))=nanmean(c(1,:));
                    end
                    if sum(isnan(c(2,:)))==numel(c(2,:))
                         c(2,find(isnan(c(2,:))))=0;
                    else
                         c(2,find(isnan(c(2,:))))=nanmean(c(2,:));
                    end

                    matlabbatch{1}.spm.stats.factorial_design.cov(1).c = c(:);
                    matlabbatch{1}.spm.stats.factorial_design.cov(1).cname = ... 
                        GLMFlags.Covariates{GLMFlags.ActiveCov(cx)};
                end

  
                 cd(matlabbatch{1}.spm.stats.factorial_design.dir{1}) ; save Design matlabbatch; cd(cwd);
                spm_jobman('run',matlabbatch)
            end
                
            for ana = 1:numel(GLMFlags.Label)         
                load(fullfile(cwd,'matfiles','se_Contrasts.mat'))
                matlabbatch{1}.spm.stats.con.spmmat{1} = ... 
                    fullfile(cwd,'Projects','RS','GLM',mgroup,GLMFlags.project,... 
                    des,folderNow,GLMFlags.Label{ana},'SPM.mat');

                matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'T1-T2';
                matlabbatch{1}.spm.stats.con.consess{1}.tcon.convec = [1 -1];

                matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'T2-T1';
                matlabbatch{1}.spm.stats.con.consess{2}.tcon.convec = [-1 1];
                
                prs=numel(GLMFlags.sub)/2;
                matlabbatch{1}.spm.stats.con.consess{3}.tcon.name = 'main effect T1+';
                matlabbatch{1}.spm.stats.con.consess{3}.tcon.convec = [1 0 0 0 ones(1,prs)/prs];
              
                matlabbatch{1}.spm.stats.con.consess{4}.tcon.name = 'main effect T1-';
                matlabbatch{1}.spm.stats.con.consess{4}.tcon.convec = [-1 0 0 0 -ones(1,prs)/prs];

                matlabbatch{1}.spm.stats.con.consess{5}.tcon.name = 'main effect T2+';
                matlabbatch{1}.spm.stats.con.consess{5}.tcon.convec = [0 1 0 0 ones(1,prs)/prs];

                matlabbatch{1}.spm.stats.con.consess{6}.tcon.name = 'main effect T2-';
                matlabbatch{1}.spm.stats.con.consess{6}.tcon.convec = [0 -1 0 0 -ones(1,prs)/prs];

                cnt=7;
                for cx = 1:numel(GLMFlags.ActiveCov)
                    matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.name = ... 
                        ['T1-T2 ' GLMFlags.Covariates{GLMFlags.ActiveCov(cx)}];
                    matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.convec = [0 0 1 -1];
                    cnt=cnt+1;
                    matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.name = ... 
                        ['T2-T1 ' GLMFlags.Covariates{GLMFlags.ActiveCov(cx)}];
                    matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.convec = [0 0 -1 1];
                    cnt=cnt+1;
                    matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.name = ... 
                        ['main effect T1: ' GLMFlags.Covariates{GLMFlags.ActiveCov(cx)} ' +'];
                    matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.convec = [0 0 1 0 ones(1,prs)/prs];
                    cnt=cnt+1;
                    matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.name = ... 
                        ['main effect T1: ' GLMFlags.Covariates{GLMFlags.ActiveCov(cx)} ' -'];
                    matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.convec = [0 0 -1 0 -ones(1,prs)/prs];
                    cnt=cnt+1;
                    matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.name = ... 
                        ['main effect T2: ' GLMFlags.Covariates{GLMFlags.ActiveCov(cx)} ' +'];
                    matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.convec = [0 0 0 1 ones(1,prs)/prs];
                    cnt=cnt+1;
                    matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.name = ... 
                        ['main effect T2: ' GLMFlags.Covariates{GLMFlags.ActiveCov(cx)} ' -'];
                    matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.convec = [0 0 0 -1 -ones(1,prs)/prs];
                    cnt=cnt+1;
                end
                
                spm_jobman('run',matlabbatch)
            end                 
    end
    
    
    
    
    
    
else
    analyses = {'MainEffect','Covariate','MultipleCovariates'};
    
    folderNow = ['ANOVA_' analyses{GLMFlags.gotAnalysis}];
    
    
    
    
    if GLMFlags.gotBlock>1
        GLMFlags.xGroup = GLMFlags.Cov(:,GLMFlags.gotBlock-1);
        if ~all(unique(GLMFlags.xGroup)==[1:numel(unique(GLMFlags.xGroup))]')
            tmp = unique(GLMFlags.xGroup);
            vals = unique(tmp);
            GLMFlags.xGroup = ones(1,numel(GLMFlags.xGroup));
            for i=1:numel(vals)
                GLMFlags.xGroup(GLMFlags.xGroup==vals(i)) == i;
            end
        end
    else
        GLMFlags.xGroup = ones(size(GLMFlags.Cov,1),1);
    end
    
    
    switch GLMFlags.gotAnalysis
        case 1
            load(fullfile(cwd,'matfiles','se_ANOVA.mat'))
            matlabbatch{1}.spm.stats.factorial_design.masking.em{1} = fullfile(cwd,'MaskenEtc','MNI152GM.nii');
            matlabbatch{1}.spm.stats.factorial_design.masking.im = 0;
            matlabbatch{1}.spm.stats.factorial_design.dir{1} = ... 
                    fullfile(cwd,'Projects','RS','GLM',mgroup,GLMFlags.project,des,folderNow);
            mkdir(matlabbatch{1}.spm.stats.factorial_design.dir{1});
            for s=1:numel(GLMFlags.sub)
                matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(s).scans = {};
                for ana = 1:numel(GLMFlags.Label)
                    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(s).scans{ana,1} = ...
                        fullfile(cwd,'Projects','RS','GLM',mgroup,GLMFlags.project, ... 
                        des,GLMFlags.Label{ana},[GLMFlags.sub{s}  '.nii']);
                end
                matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(s).conds = ... 
                    [1:numel(GLMFlags.Label)];
            end
            %      matlabbatch = se_adjustAnova(matlabbatch);
             cd(matlabbatch{1}.spm.stats.factorial_design.dir{1}) ; save Design matlabbatch; cd(cwd);
             spm_jobman('run',matlabbatch)
            
            
            load(fullfile(cwd,'matfiles','se_Contrasts.mat'))
            matlabbatch{1}.spm.stats.con.consess = {matlabbatch{1}.spm.stats.con.consess{1}};
            matlabbatch{1}.spm.stats.con.spmmat{1} = ... 
                fullfile(cwd,'Projects','RS','GLM',mgroup,GLMFlags.project,des,folderNow,'SPM.mat');
            matlabbatch{1}.spm.stats.con.consess{1}.fcon.name = 'EoI';
            try; matlabbatch{1}.spm.stats.con.consess{1} = ... 
                rmfield(matlabbatch{1}.spm.stats.con.consess{1},'tcon'); end
            matlabbatch{1}.spm.stats.con.consess{1}.fcon.convec = {eye(numel(GLMFlags.Label))};
            spm_jobman('run',matlabbatch)
            
            load(fullfile(cwd,'matfiles','se_Contrasts.mat'))
            matlabbatch{1}.spm.stats.con.spmmat{1} = ... 
                fullfile(cwd,'Projects','RS','GLM',mgroup,GLMFlags.project,des,folderNow,'SPM.mat');
            for ana = 1:numel(GLMFlags.Label)
                matlabbatch{1}.spm.stats.con.consess{ana}.tcon.name = [GLMFlags.Label{ana} ' +'];
                vec = zeros(1,numel(GLMFlags.Label)); vec(ana)  = 1;
                matlabbatch{1}.spm.stats.con.consess{ana}.tcon.convec = vec;
            end
            for ana = 1:numel(GLMFlags.Label)
                matlabbatch{1}.spm.stats.con.consess{ana+(1*numel(GLMFlags.Label))}.tcon.name = ... 
                    [GLMFlags.Label{ana} ' -'];
                vec = zeros(1,numel(GLMFlags.Label)); vec(ana)  = -1;
                matlabbatch{1}.spm.stats.con.consess{ana+(1*numel(GLMFlags.Label))}.tcon.convec = vec;
            end
            spm_jobman('run',matlabbatch)
            
            
        case 2
            xfolder = folderNow;
            
            for cx = 1:numel(GLMFlags.ActiveCov)
                load(fullfile(cwd,'matfiles','se_ANOVAcov.mat'))
                matlabbatch{1}.spm.stats.factorial_design.masking.em{1} = fullfile(cwd,'MaskenEtc','MNI152GM.nii');
                matlabbatch{1}.spm.stats.factorial_design.masking.im = 0;
                folderNow = [xfolder '_' GLMFlags.Covariates{GLMFlags.ActiveCov(cx)}];
                matlabbatch{1}.spm.stats.factorial_design.dir{1} = ... 
                    fullfile(cwd,'Projects','RS','GLM',mgroup,GLMFlags.project,des,folderNow);
                mkdir(matlabbatch{1}.spm.stats.factorial_design.dir{1});
                for s=1:numel(GLMFlags.sub)
                    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(s).scans = {};
                    for ana = 1:numel(GLMFlags.Label)
                        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(s).scans{ana,1} = ...
                            fullfile(cwd,'Projects','RS','GLM',mgroup,GLMFlags.project, ... 
                            des,GLMFlags.Label{ana},[GLMFlags.sub{s}  '.nii']);
                    end
                    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(s).conds = ... 
                        [1:numel(GLMFlags.Label)];
                end
                matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject = ...
                    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject ... 
                    (find(sum(isnan(GLMFlags.Cov(:,GLMFlags.ActiveCov(cx))),2)==0));
                
                c = GLMFlags.Cov(:,GLMFlags.ActiveCov(cx));
                c =c(find(sum(isnan(GLMFlags.Cov(:,GLMFlags.ActiveCov(cx))),2)==0));
                c = repmat(c',numel(GLMFlags.Label),1);
                matlabbatch{1}.spm.stats.factorial_design.cov(1).c = c (:);
                matlabbatch{1}.spm.stats.factorial_design.cov(1).cname = ... 
                    GLMFlags.Covariates{GLMFlags.ActiveCov(cx)};
                cd(matlabbatch{1}.spm.stats.factorial_design.dir{1}) ; save Design matlabbatch; cd(cwd);
                spm_jobman('run',matlabbatch)
                
                
                load(fullfile(cwd,'matfiles','se_Contrasts.mat'))
                matlabbatch{1}.spm.stats.con.spmmat{1} = ... 
                    fullfile(cwd,'Projects','RS','GLM',mgroup,GLMFlags.project,des,folderNow,'SPM.mat');
                matlabbatch{1}.spm.stats.con.consess{1}.fcon.name = 'EoI';
                try; matlabbatch{1}.spm.stats.con.consess{1} = ... 
                    rmfield(matlabbatch{1}.spm.stats.con.consess{1},'tcon'); end
                matlabbatch{1}.spm.stats.con.consess{1}.fcon.convec = {eye(numel(GLMFlags.Label)*2)};
                matlabbatch{1}.spm.stats.con.consess{2}.fcon.name = 'EoConnectivity';
                try; matlabbatch{1}.spm.stats.con.consess{2} = ... 
                    rmfield(matlabbatch{1}.spm.stats.con.consess{2},'tcon'); end
                matlabbatch{1}.spm.stats.con.consess{2}.fcon.convec = ... 
                    {[eye(numel(GLMFlags.Label)) zeros(numel(GLMFlags.Label))]};
                matlabbatch{1}.spm.stats.con.consess{3}.fcon.name = 'EoCovariate';
                try; matlabbatch{1}.spm.stats.con.consess{2} = ... 
                    rmfield(matlabbatch{1}.spm.stats.con.consess{3},'tcon'); end
                matlabbatch{1}.spm.stats.con.consess{3}.fcon.convec = ... 
                    {[zeros(numel(GLMFlags.Label)) eye(numel(GLMFlags.Label))]};
                spm_jobman('run',matlabbatch)
                
                load(fullfile(cwd,'matfiles','se_Contrasts.mat'))
                matlabbatch{1}.spm.stats.con.spmmat{1} = ... 
                    fullfile(cwd,'Projects','RS','GLM',mgroup,GLMFlags.project,des,folderNow,'SPM.mat');
                for ana = 1:numel(GLMFlags.Label)
                    matlabbatch{1}.spm.stats.con.consess{ana}.tcon.name = [GLMFlags.Label{ana} ' +'];
                    vec = zeros(1,numel(GLMFlags.Label)*2); vec(ana)  = 1;
                    matlabbatch{1}.spm.stats.con.consess{ana}.tcon.convec = vec;
                end
                for ana = 1:numel(GLMFlags.Label)
                    matlabbatch{1}.spm.stats.con.consess{ana+(1*numel(GLMFlags.Label))}.tcon.name = ... 
                        [GLMFlags.Label{ana} ' -'];
                    vec = zeros(1,numel(GLMFlags.Label)*2); vec(ana)  = -1;
                    matlabbatch{1}.spm.stats.con.consess{ana+(1*numel(GLMFlags.Label))}.tcon.convec = vec;
                end
                for ana = 1:numel(GLMFlags.Label)
                    matlabbatch{1}.spm.stats.con.consess{ana+(2*numel(GLMFlags.Label))}.tcon.name = ... 
                        [GLMFlags.Label{ana} ': ' GLMFlags.Covariates{GLMFlags.ActiveCov(cx)} ' +'];
                    vec = zeros(1,numel(GLMFlags.Label)*2); vec(ana+numel(GLMFlags.Label))  = 1;
                    matlabbatch{1}.spm.stats.con.consess{ana+(2*numel(GLMFlags.Label))}.tcon.convec = vec;
                end
                for ana = 1:numel(GLMFlags.Label)
                    matlabbatch{1}.spm.stats.con.consess{ana+(3*numel(GLMFlags.Label))}.tcon.name = ... 
                        [GLMFlags.Label{ana} ': ' GLMFlags.Covariates{GLMFlags.ActiveCov(cx)} ' -'];
                    vec = zeros(1,numel(GLMFlags.Label)*2); vec(ana+numel(GLMFlags.Label))  = -1;
                    matlabbatch{1}.spm.stats.con.consess{ana+(3*numel(GLMFlags.Label))}.tcon.convec = vec;
                end
                spm_jobman('run',matlabbatch)
                
            end
            
        case 3
            load(fullfile(cwd,'matfiles','se_ANOVAcov.mat'))
            matlabbatch{1}.spm.stats.factorial_design.masking.em{1} = fullfile(cwd,'MaskenEtc','MNI152GM.nii');
            matlabbatch{1}.spm.stats.factorial_design.masking.im = 0;
            matlabbatch{1}.spm.stats.factorial_design.dir{1} = ... 
                fullfile(cwd,'Projects','RS','GLM',mgroup,GLMFlags.project,des,folderNow);
            mkdir(matlabbatch{1}.spm.stats.factorial_design.dir{1});
            for s=1:numel(GLMFlags.sub)
                matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(s).scans = {};
                for ana = 1:numel(GLMFlags.Label)
                    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(s).scans{ana,1} = ...
                        fullfile(cwd,'Projects','RS','GLM',mgroup,GLMFlags.project,... 
                        des,GLMFlags.Label{ana},[GLMFlags.sub{s}  '.nii']);
                end
                matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(s).conds = ... 
                    [1:numel(GLMFlags.Label)];
            end
            matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject =...
                matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject... 
                    (find(sum(isnan(GLMFlags.Cov(:,GLMFlags.ActiveCov)),2)==0));
            
            for cx = 1:numel(GLMFlags.ActiveCov)
                matlabbatch{1}.spm.stats.factorial_design.cov(cx) = ... 
                    matlabbatch{1}.spm.stats.factorial_design.cov(1);
                c = GLMFlags.Cov(:,GLMFlags.ActiveCov(cx));
                c =c(find(sum(isnan(GLMFlags.Cov(:,GLMFlags.ActiveCov)),2)==0));
                c = repmat(c',numel(GLMFlags.Label),1);
                matlabbatch{1}.spm.stats.factorial_design.cov(cx).c = c (:);
                matlabbatch{1}.spm.stats.factorial_design.cov(cx).cname = ... 
                    GLMFlags.Covariates{GLMFlags.ActiveCov(cx)};
            end
            %      matlabbatch = se_adjustAnova(matlabbatch);
             cd(matlabbatch{1}.spm.stats.factorial_design.dir{1}) ; save Design matlabbatch; cd(cwd);
             spm_jobman('run',matlabbatch)
            
            
            load(fullfile(cwd,'matfiles','se_Contrasts.mat')); cnt = 1;
            matlabbatch{1}.spm.stats.con.spmmat{1} = ... 
                fullfile(cwd,'Projects','RS','GLM',mgroup,GLMFlags.project,des,folderNow,'SPM.mat');
            for ana = 1:numel(GLMFlags.Label)
                matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.name = [GLMFlags.Label{ana} ' +'];
                vec = zeros(1,numel(GLMFlags.Label)*(numel(GLMFlags.ActiveCov)+1));
                vec(ana)  = 1;  matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.convec = vec; cnt = cnt+1;
                matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.name = [GLMFlags.Label{ana} ' -'];
                vec = zeros(1,numel(GLMFlags.Label)*(numel(GLMFlags.ActiveCov)+1));
                vec(ana)  = -1; matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.convec = vec; cnt = cnt+1;
            end
            spm_jobman('run',matlabbatch)
            
            load(fullfile(cwd,'matfiles','se_Contrasts.mat'))
            matlabbatch{1}.spm.stats.con.spmmat{1} = ... 
                fullfile(cwd,'Projects','RS','GLM',mgroup,GLMFlags.project,des,folderNow,'SPM.mat');
            cnt = 1;
            for ana = 1:numel(GLMFlags.Label)
                for i=1:numel(GLMFlags.ActiveCov)
                    matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.name = ... 
                        [GLMFlags.Label{ana} ': ' GLMFlags.Covariates{GLMFlags.ActiveCov(i)} ' +'];
                    vec = zeros(1,numel(GLMFlags.Label)*(numel(GLMFlags.ActiveCov)+1));
                    vec(((i)*numel(GLMFlags.Label))+ana)  = 1;
                    matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.convec = vec;
                    cnt = cnt+1;
                    
                    matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.name = ... 
                    [GLMFlags.Label{ana} ': ' GLMFlags.Covariates{GLMFlags.ActiveCov(i)}  ' -'];
                    vec = zeros(1,numel(GLMFlags.Label)*(numel(GLMFlags.ActiveCov)+1));
                    vec(((i)*numel(GLMFlags.Label))+ana)  = -1;
                    matlabbatch{1}.spm.stats.con.consess{cnt}.tcon.convec = vec;
                    cnt = cnt+1;
                end
            end
            spm_jobman('run',matlabbatch)
            
            load(fullfile(cwd,'matfiles','se_Contrasts.mat'))
            matlabbatch{1}.spm.stats.con.spmmat{1} = ... 
                fullfile(cwd,'Projects','RS','GLM',mgroup,GLMFlags.project,des,folderNow,'SPM.mat');
            matlabbatch{1}.spm.stats.con.consess{1}.fcon.name = 'EoConnectivity';
            try; matlabbatch{1}.spm.stats.con.consess{1} = ... 
                rmfield(matlabbatch{1}.spm.stats.con.consess{1},'tcon'); end
            matlabbatch{1}.spm.stats.con.consess{1}.fcon.convec = ... 
                {[eye(numel(GLMFlags.Label)) repmat(zeros(numel(GLMFlags.Label)),1,numel(GLMFlags.ActiveCov))]};
            matlabbatch{1}.spm.stats.con.consess{2}.fcon.name = 'EoCovariates';
            try; matlabbatch{1}.spm.stats.con.consess{2} = ... 
                rmfield(matlabbatch{1}.spm.stats.con.consess{2},'tcon'); end
            matlabbatch{1}.spm.stats.con.consess{2}.fcon.convec = ... 
                {[zeros(numel(GLMFlags.Label)*numel(GLMFlags.ActiveCov),numel(GLMFlags.Label)) ... 
                    eye(numel(GLMFlags.Label)*numel(GLMFlags.ActiveCov))]};
            
            cnt = 3;
            for cx=1:numel(GLMFlags.ActiveCov)
                matlabbatch{1}.spm.stats.con.consess{cnt}.fcon.name = ... 
                    ['Eo' GLMFlags.Covariates{GLMFlags.ActiveCov(cx)}];
                try; matlabbatch{1}.spm.stats.con.consess{cnt} = ... 
                    rmfield(matlabbatch{1}.spm.stats.con.consess{cnt},'tcon'); end
                vec = zeros(numel(GLMFlags.Label),numel(GLMFlags.Label)*(numel(GLMFlags.ActiveCov)+1));
                for ana = 1:numel(GLMFlags.Label); vec(ana,((cx)*numel(GLMFlags.Label))+ana)  = 1; end
                matlabbatch{1}.spm.stats.con.consess{cnt}.fcon.convec = {vec};
                cnt = cnt+1;
            end
            spm_jobman('run',matlabbatch)
            
            
    end
    
end


for ana = 1:numel(GLMFlags.Label)
    
    zip(fullfile(cwd,'Projects','RS','GLM',mgroup,GLMFlags.project,des,[GLMFlags.Label{ana} '.zip']),... 
            fullfile(cwd,'Projects','RS','GLM',mgroup,GLMFlags.project,des,GLMFlags.Label{ana}))
    rmdir(fullfile(cwd,'Projects','RS','GLM',mgroup,GLMFlags.project,des,GLMFlags.Label{ana}),'s')
end

