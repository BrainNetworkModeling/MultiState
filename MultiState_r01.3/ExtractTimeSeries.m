clc, clear

fid = fopen('_path2data','r'); datadir = fgetl(fid); fclose(fid);

NetFlags.Group = spm_select(1,'any',['Select group description'],[],fullfile(pwd,'Lookups'),'.mat',1);
tmp = load(NetFlags.Group);
NetFlags.Covariates = tmp.Covariates;
NetFlags.sub = tmp.sub(tmp.hasRS==1 & tmp.hasVBM==1);
NetFlags.Cov = tmp.Cov(tmp.hasRS==1 & tmp.hasVBM==1,:);
NetFlags.isPat = tmp.isPat(tmp.hasRS==1 & tmp.hasVBM==1);
NetFlags.SubDir = tmp.SubDir(tmp.hasRS==1 & tmp.hasVBM==1);
clear tmp


TxtName  = spm_select(1,'any',['Select coordinate or VOI list'],[],fullfile(pwd,'VOIs'),'.txt',1);
try
    [NetFlags.X NetFlags.Y NetFlags.Z NetFlags.Label]     = textread(TxtName,'%f %f %f %s');
    NetFlags.Coordinates = 1;
catch
    [NetFlags.VOI NetFlags.Label] = textread(TxtName,'%s %s');
    for i=1:numel(NetFlags.VOI)
        try
            spm_vol(NetFlags.VOI{1});
        catch
            spm_vol(fullfile(pwd,'VOIs','VOIfiles',spm_str_manip(NetFlags.VOI{1},'t')));
            NetFlags.VOI{1} = fullfile(pwd,'VOIs','VOIfiles',spm_str_manip(NetFlags.VOI{1},'t'));
        end
    end
    NetFlags.Coordinates = 0;
end
NetFlags.project = strrep(strrep(spm_str_manip(TxtName,'rt'),'_VOIs',''),'_ImageList','');



time = zeros(1,numel(NetFlags.sub));
TR = zeros(1,numel(NetFlags.sub)); 
scans = zeros(1,numel(NetFlags.sub)); 
for s=1:numel(NetFlags.sub)
    tmp = load(fullfile(datadir,'DATA',NetFlags.SubDir{s},NetFlags.sub{s},'Confounds5mmPCA.mat'));
    time(s) = numel(tmp.V)*tmp.TR/60;
    TR(s) = tmp.TR;
    scans(s) = numel(tmp.V);
end
clear tmp; cut = prctile(time,10);

NetFlags.TR    = mode(TR(time>=cut));
NetFlags.scans = mode(scans(find(time==min(time(time>=cut)))));
NetFlags.time =  mode(time(find(time==min(time(time>=cut)))));


NetFlags.sub = NetFlags.sub(time >= cut);
NetFlags.Cov = NetFlags.Cov(time >= cut,:);
NetFlags.isPat = NetFlags.isPat(time >= cut);
NetFlags.SubDir = NetFlags.SubDir(time >= cut);




NetFlags.Mask = 1;
NetFlags.PCA = 2;
NetFlags.Glob = 4;
NetFlags.Radius = 5;


GM   = {'NoGMmask','IndGMmask','GroupGMmask'};
PCA  = {'PCA','NoPCA'};
Glob = {'TSR','lTSR','WMCSF','lWMCSF','GSR','lGSR','NoGSR'};
if NetFlags.Coordinates==0
    des = [GM{NetFlags.Mask} '_' PCA{NetFlags.PCA} '_' Glob{NetFlags.Glob}];
else
    des = [GM{NetFlags.Mask} '_' PCA{NetFlags.PCA}  '_' Glob{NetFlags.Glob} '_' strrep(num2str(NetFlags.Radius,'%2.1f'),'.','')];
end

            
Qreg = 9:32;
if NetFlags.PCA==1; Qreg = [Qreg 33:37]; end
if NetFlags.Glob==1;      Qreg = [Qreg 1:8];
elseif NetFlags.Glob==2;  Qreg = [Qreg 1:4];
elseif NetFlags.Glob==3;  Qreg = [Qreg 2 3 6 7];
elseif NetFlags.Glob==4;  Qreg = [Qreg 2:3];    
elseif NetFlags.Glob==5;  Qreg = [Qreg 4 8];
elseif NetFlags.Glob==6;  Qreg = [Qreg 4];
end

[~, ~, ~] = mkdir(fullfile(pwd,'Projects','RS','ROIextract',spm_str_manip(NetFlags.Group,'rt')));
load Maske


 allData = single(nan(numel(NetFlags.sub),NetFlags.scans,numel(NetFlags.Label)));
 allVol  = single(nan(numel(NetFlags.sub),numel(NetFlags.Label)));


VOI  = [];
Ind  = [];

for ana = 1:numel(NetFlags.Label)
    if NetFlags.Coordinates==0
        try
            Vx = spm_vol(NetFlags.VOI{ana});
        catch
            try
                Vx = spm_vol(fullfile(pwd,'VOIs','VOIfiles',spm_str_manip(NetFlags.VOI{ana},'t')));
            catch
                Vx = spm_vol(fullfile(pwd,'VOIfiles',spm_str_manip(NetFlags.VOI{ana},'t')));
            end
        end
        sXYZ = Vx.mat \ VM.mat * [maskXYZ; ones(1,size(maskXYZ,2))];
        D = spm_get_data(Vx,sXYZ);
        Ind = [Ind find(D>0)];
        VOI = [VOI ones(1,sum(D>0))*ana];
    
        dat = spm_read_vols(Vx);
        tmp = find(dat); Zdat{ana} = dat(tmp); clear XYZ
        [XYZ(:,1) XYZ(:,2) XYZ(:,3)] = ind2sub(Vx.dim,tmp);
        tXYZ{ana} = (Vx.mat) * [XYZ'; ones(1,size(XYZ,1))];
        allVOIs(ana) = Vx;

    else
        sXYZ = VM.mat * [maskXYZ; ones(1,size(maskXYZ,2))];
        D    = sqrt((sXYZ(1,:)-NetFlags.X(ana)).^2+(sXYZ(2,:)-NetFlags.Y(ana)).^2+(sXYZ(3,:)-NetFlags.Z(ana)).^2);
        Ind = [Ind find(D<NetFlags.Radius)];
        VOI = [VOI ones(1,sum(D<NetFlags.Radius))*ana];
    end
    
end

if NetFlags.Mask==3
    Vgm = spm_vol(fullfile(pwd,'MaskenEtc','forRS',[spm_str_manip(NetFlags.Group,'rt') '.nii']));
    gm  = spm_get_data(Vgm,Vgm.mat \ VM.mat * [maskXYZ(:,Ind); ones(1,size(Ind,2))]);
elseif NetFlags.Mask==1
    gm = ones(1,numel(Ind));
end


for s=1:numel(NetFlags.sub)
    fprintf(1,'%s\n',['Subject ' int2str(s) ' / ' int2str(numel(NetFlags.sub))]);
    
    load(fullfile(datadir,'DATA',NetFlags.SubDir{s},NetFlags.sub{s},'Confounds5mmPCA.mat'));
    for i=1:numel(V); V(i).fname = fullfile(datadir,'DATA',NetFlags.SubDir{s},NetFlags.sub{s},'RS',V(i).fname); end
    regnow = [reg(:,Qreg) ones(size(reg,1),1)];
    
    if NetFlags.Mask==2
        try
            Vgm = spm_vol(fullfile(datadir,'DATA',NetFlags.SubDir{s},NetFlags.sub{s},['c1w' NetFlags.sub{s} '_meanEPI.nii']));
        catch
            Vgm = spm_vol(fullfile(pwd,'MaskenEtc','forRS',[spm_str_manip(NetFlags.Group,'rt') '.nii']));
        end
        gm  = spm_get_data(Vgm,Vgm.mat \ VM.mat * [maskXYZ(:,Ind); ones(1,size(Ind,2))]);
    end
    
    QQ = [];
    for ana = 1:numel(NetFlags.Label)
        xQ = find(VOI==ana);
        QQ = [QQ xQ(gm(xQ)>=median(gm(xQ)))];
    end
    
    XYZ = [Ind(QQ); ones(2,numel(QQ))];
    y    = zeros(numel(V),size(XYZ,2));
    for i = 1:length(V)
        y(i,:) = spm_sample_vol(V(i),XYZ(1,:),XYZ(2,:),XYZ(3,:),0);
    end
    
    bb = regnow\y;
    y = y-regnow*bb; clear bb
    y = conn_filter(TR,[0.01, 0.08],y);
    
    clear TS
    
    for ana = 1:numel(NetFlags.Label)
        try
            TS(:,ana)  = ev(y(:,VOI(QQ)==ana));
        catch
            TS(:,ana)  = mean(y(:,VOI(QQ)==ana),2);
            TS(:,ana)  = TS(:,ana) - mean(TS(:,ana));
            fprintf(1,'%s\n',['*********         *********     Subject ' int2str(s)]);
        end
    end
    
    TS = single(interp1(TR:TR:numel(V)*TR,TS,NetFlags.TR:NetFlags.TR:NetFlags.time*60));
    allData(s,:,:) = TS;
    
    
    fil = dir(fullfile(datadir,'DATA',NetFlags.SubDir{s},NetFlags.sub{s},'3D',[ 'm0wrp1*.nii']));
    V1 = spm_vol(fullfile(datadir,'DATA',NetFlags.SubDir{s},NetFlags.sub{s},'3D',fil(1).name));
    
    
    for ana = 1:numel(NetFlags.Label)
        allVol(s,ana) = nansum(spm_get_data(V1,inv(V1.mat) * tXYZ{ana}).*Zdat{ana}')*abs(det(allVOIs(ana).mat));
    end
    
end

allData = allData(:,2:end-1,:);



NetFlags = rmfield(NetFlags,'sub');
NetFlags = rmfield(NetFlags,'SubDir');
NetFlags = rmfield(NetFlags,'Cov');
NetFlags = rmfield(NetFlags,'Covariates');

save(fullfile(pwd,'Projects','RS','ROIextract',spm_str_manip(NetFlags.Group,'rt'),[NetFlags.project '.mat']),'allData','allVol','NetFlags')