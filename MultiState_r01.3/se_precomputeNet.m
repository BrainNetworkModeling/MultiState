function se_precomputeNet(des)

global NetFlags
global fg
cwd=pwd;

fid = fopen('_path2data','r'); datadir = fgetl(fid); fclose(fid);

try; NetFlags.gotMask = get(NetFlags.Mask,'value'); end
try; NetFlags.gotPCA = get(NetFlags.PCA,'value'); end
try; NetFlags.gotGlob = get(NetFlags.Glob,'value'); end
try; NetFlags.gotFilt = get(NetFlags.Filt,'value'); end

    GM   = {'NoGMmask','IndGMmask','GroupGMmask'};
    PCA  = {'PCA','NoPCA','FIX'};
    Glob = {'TSR','lTSR','WMCSF','lWMCSF','GSR','lGSR','NoGSR'};
    FILT = {'BP','HP'};

if NetFlags.Coordinates==0
    des = [GM{NetFlags.gotMask} '_' PCA{NetFlags.gotPCA} '_' Glob{NetFlags.gotGlob}...
        '_' FILT{NetFlags.gotFilt}];
else
    des = [GM{NetFlags.gotMask} '_' PCA{NetFlags.gotPCA}  '_' Glob{NetFlags.gotGlob}...
        '_' FILT{NetFlags.gotFilt} '_' strrep(num2str(NetFlags.Radius,'%2.1f'),'.','')];
end

if NetFlags.gotPCA==3
    RSdir='RS_fix';
    Qreg = [];       % 24 movement regressors already corrected for in the FIX denoising
else
    RSdir='RS';
    Qreg = 9:32;    % 24 movement regressors from the realignment parameter
end

%        Nuisance regressors as estimated in the preprocessing:
%        Nr     4    4       6  6      6    6        5
%        No     1-   5-      9- 15-   21-  27-      33-      37
%        reg = [gx2' gx2'.^2 rp rp.^2 rpp rpp.^2 COEFF(:,1:5)];

if NetFlags.gotPCA==1
    Qreg = [Qreg 33:37];
end
if NetFlags.gotGlob==1
    Qreg = [Qreg 1 2 3 5 6 7];
elseif NetFlags.gotGlob==2
    Qreg = [Qreg 1:3];
elseif NetFlags.gotGlob==3
    Qreg = [Qreg 2 3 6 7];
elseif NetFlags.gotGlob==4
    Qreg = [Qreg 2:3];
elseif NetFlags.gotGlob==5
    Qreg = [Qreg 4 8];
elseif NetFlags.gotGlob==6
    Qreg = [Qreg 4];
end

if NetFlags.gotFilt==1
    Filt = [0.01, 0.08];
else
    Filt = [0.01, inf];
end    
    
[~, ~, ~] = mkdir(fullfile(cwd,'Projects','RS','Networks',spm_str_manip(NetFlags.Group,'rt') ...
    ,[NetFlags.project '_' des]));
fid = fopen(fullfile(cwd,'Projects','RS','Networks',spm_str_manip(NetFlags.Group,'rt') ...
    ,[NetFlags.project '_' des],[NetFlags.project '_' des '_subjects.txt']),'w+');
fprintf(fid,'%s\t','Subject');

for i=1:numel(NetFlags.Covariates)
    fprintf(fid,'%s\t',NetFlags.Covariates{i});
end
fprintf(fid,'%s\t','is Patient');
for ii=1:numel(NetFlags.sub)
    fprintf(fid,'\n%s\t',NetFlags.sub{ii});
    
    for i=1:numel(NetFlags.Covariates)
        fprintf(fid,'%4.1f\t',NetFlags.Cov(ii,i));
    end
        fprintf(fid,'%4.0f\t',NetFlags.isPat(ii));
end
fclose(fid);


load ./MaskenEtc/MNI152GM.mat
NN = [];

z = nan(numel(NetFlags.sub),numel(NetFlags.Label),numel(NetFlags.Label));

fullfile(spm_str_manip(NetFlags.Group,'rt'),[NetFlags.project '_' des],[NetFlags.project '_' des '.mat'])


VOI  = [];
Ind  = [];

for ana = 1:numel(NetFlags.Label)
    
    if NetFlags.Coordinates==0
        try
            Vx = spm_vol(NetFlags.VOI{ana});
        catch
            try
                Vx = spm_vol(fullfile(cwd,'VOIs','VOIfiles',spm_str_manip(NetFlags.VOI{ana},'t')));
            catch
                Vx = spm_vol(fullfile(cwd,'VOIfiles',spm_str_manip(NetFlags.VOI{ana},'t')));
            end
        end
        sXYZ = Vx.mat \ VM.mat * [maskXYZ; ones(1,size(maskXYZ,2))];
        D = spm_get_data(Vx,sXYZ);
        Ind = [Ind find(D>0)];
        VOI = [VOI ones(1,sum(D>0))*ana];
    else
        sXYZ = VM.mat * [maskXYZ; ones(1,size(maskXYZ,2))];
        D   = sqrt((sXYZ(1,:)-NetFlags.X(ana)).^2+(sXYZ(2,:)-NetFlags.Y(ana)).^2+(sXYZ(3,:)-NetFlags.Z(ana)).^2);
        Ind = [Ind find(D<NetFlags.Radius)];
        VOI = [VOI ones(1,sum(D<NetFlags.Radius))*ana];
    end
    
end

if NetFlags.gotMask==3
    Vgm = spm_vol(fullfile(cwd,'MaskenEtc','forRS',[spm_str_manip(NetFlags.Group,'rt') '.nii']));
    gm  = spm_get_data(Vgm,Vgm.mat \ VM.mat * [maskXYZ(:,Ind); ones(1,size(Ind,2))]);
elseif NetFlags.gotMask==1
    gm = ones(1,numel(Ind));
end

[~, ~, ~] = mkdir(fullfile(cwd,'Projects','RS','Networks',spm_str_manip(NetFlags.Group,'rt'), ...
[NetFlags.project '_' des],'TS','TS_all'));    % folder for TS

if any(NetFlags.isPat)>0
    pat=0;con=0;
    [~, ~, ~] = mkdir(fullfile(cwd,'Projects','RS','Networks',spm_str_manip(NetFlags.Group,'rt') ...
        ,[NetFlags.project '_' des],'TS','TS_con'));    % folder for TS
    [~, ~, ~] = mkdir(fullfile(cwd,'Projects','RS','Networks',spm_str_manip(NetFlags.Group,'rt') ...
        ,[NetFlags.project '_' des],'TS','TS_pat'));    % folder for TS
end

if 1==1 % more than 1 site
    fprintf(1,'%s\n',['checking number of timepoints for ' int2str(numel(NetFlags.sub)) ' subjects']);
    for s=1:numel(NetFlags.sub)
        tmp = load(fullfile(datadir,'DATA',NetFlags.SubDir{s},NetFlags.sub{s},RSdir, ...
            ['Counfounds_' NetFlags.sub{s} '.mat']));
        time(s) = numel(tmp.Vi)*tmp.TR/60;
        TR(s) = tmp.TR;
        scans(s) = numel(tmp.Vi);
        if rem(s,10)==0 
            fprintf(1,'%s','.');
        end
    end

    clear tmp; cut = prctile(time,10);

    tr    = mode(TR(time>=cut));
    NetFlags.scans = mode(scans(find(time==min(time(time>=cut)))));
    NetFlags.time =  mode(time(find(time==min(time(time>=cut)))));
    
    fprintf(1,'\n\n%s\n',['majority of time points: ' int2str(NetFlags.scans)]);
    fprintf(1,'%s\n',['maximal time points: ' int2str(max(scans))]);
    fprintf(1,'%s\n',['minimal time points: ' int2str(min(scans))]);


%     if any(scans<NetFlags.scans)    
%         fprintf(1,'\n%s\n',[num2str(numel(find(scans<NetFlags.scans))) ' subjects with less then ' ...
%             int2str(NetFlags.scans) ' time points are excluded']);
%         NetFlags.sub = NetFlags.sub(scans>=NetFlags.scans);
%         NetFlags.Cov = NetFlags.Cov(scans>=NetFlags.scans);
%         NetFlags.isPat = NetFlags.isPat(scans>=NetFlags.scans);
%         NetFlags.SubDir = NetFlags.SubDir(scans>=NetFlags.scans);
%     end
    
%     NetFlags.sub = NetFlags.sub(time >= cut);
%     NetFlags.Cov = NetFlags.Cov(time >= cut,:);
%     NetFlags.isPat = NetFlags.isPat(time >= cut);
%     NetFlags.SubDir = NetFlags.SubDir(time >= cut);

end

for s=1:numel(NetFlags.sub)
    fprintf(1,'%s\n',['Subject ' int2str(s) ' / ' int2str(numel(NetFlags.sub))]);
    
    load(fullfile(datadir,'DATA',NetFlags.SubDir{s},NetFlags.sub{s},RSdir,['Counfounds_' NetFlags.sub{s} '.mat']));
    V.fname = fullfile(datadir,'DATA',NetFlags.SubDir{s},NetFlags.sub{s},RSdir,spm_str_manip(V.fname,'t'));
    regnow = [reg(:,Qreg) ones(size(reg,1),1)];
    y = spm_read_vols(V);
    
    y = y(:,Ind);
    
    bb = regnow\y;
    y = y-regnow*bb; clear bb
    y = conn_filter(TR,Filt,y);
    
    if NetFlags.gotMask==2
        try
            Vgm = spm_vol(fullfile(datadir,'DATA',NetFlags.SubDir{s},NetFlags.sub{s},RSdir, ...
                ['c1w' NetFlags.sub{s} '_meanEPI.nii']));
        catch
            Vgm = spm_vol(fullfile(cwd,'MaskenEtc','forRS',[spm_str_manip(NetFlags.Group,'rt') '.nii']));
        end
        gm  = spm_get_data(Vgm,Vgm.mat \ VM.mat * [maskXYZ(:,Ind); ones(1,size(Ind,2))]);
    end
    
    clear TS
    for ana = 1:numel(NetFlags.Label)
        xQ = find(VOI==ana);
        try
            TS(:,ana)  = ev(y(:,xQ(gm(xQ)>=median(gm(xQ)))));
        catch
            TS(:,ana)  = mean(y(:,xQ(gm(xQ)>=median(gm(xQ)))),2);
            TS(:,ana)  = TS(:,ana) - mean(TS(:,ana));
            fprintf(1,'%s\n',['*********         *********     Subject ' int2str(s)]);
        end
    end
  
    
%%%   For different TRs interpolate here!

    if TR~=tr
        fprintf(1,'%s\n',['time points: ' int2str(NetFlags.scans) ' / TR:' num2str(tr)]);
        fprintf(1,'%s\n',['divergent time points: ' int2str(numel(TS(:,1))) ' / TR:' num2str(TR)]);
        TS = single(interp1(TR:TR:numel(Vi)*TR,TS,tr:tr:NetFlags.time*60+2));
        fprintf(1,'%s\n',['interpolated time points: ' int2str(numel(TS(:,1))) ' / TR:' num2str(tr)]);
        TS = TS(2:min(scans)+1,:);
    else
        TS = TS(1:min(scans),:);
    end
	
    r     = corr(TS(:,:));
    rr    = corr(TS(:,:),'type','Spearman');
    
    z(s,:,:) = atanh(r);
    zr(s,:,:) = atanh(rr);

    
%%%   write out time series as txt-files for FSLNets

%     TS = TS(1:NetFlags.scans,:);

    if s <= 9; no = '000'; elseif s <=99;  no = '00'; elseif s <=999; no = '0'; else no = ''; end
    fid = fopen(fullfile(cwd,'Projects','RS','Networks',spm_str_manip(NetFlags.Group,'rt'), ...
        [NetFlags.project '_' des],'TS','TS_all',['sub_' no int2str(s) '_ts.txt']),'w+');
    for i=1:min(scans)
        fprintf(fid,'%1.8f\t',TS(i,:));
        fprintf(fid,'\n');
    end
    fclose(fid);
    
    if any(NetFlags.isPat)>0
 
        if NetFlags.isPat(s)<0.5
            con = con+1;    if con <= 9; no = '000'; elseif  con <= 99; no = '00';
                            elseif  con <= 999; no = '0'; else    no = ''; end
            fid2 = fopen(fullfile(cwd,'Projects','RS','Networks',spm_str_manip(NetFlags.Group,'rt') ...
                ,[NetFlags.project '_' des],'TS','TS_con',['con_' no int2str(con) '_ts.txt']),'w+');
            for i=1:min(scans)
                fprintf(fid2,'%1.8f\t',TS(i,:));
                fprintf(fid2,'\n');
            end
            fclose(fid2);
        else
            pat = pat+1;    if pat <= 9; no = '000'; elseif  pat <= 99; no = '00';
                            elseif  pat <= 999; no = '0'; else    no = ''; end
            fid3 = fopen(fullfile(cwd,'Projects','RS','Networks',spm_str_manip(NetFlags.Group,'rt') ...
                ,[NetFlags.project '_' des],'TS','TS_pat',['pat_' no int2str(pat) '_ts.txt']),'w+');
            for i=1:min(scans)
                fprintf(fid3,'%1.8f\t',TS(i,:));
                fprintf(fid3,'\n');
            end
            fclose(fid3);
        end
    end
  
end

if ~exist('tr','var')
    tr=TR;
end

z(z==Inf)=0;
zr(zr==Inf)=0;

% NetFlags.z = z;
% NetFlags.rawz = z;

%%% print analysis settings
%print(fg,'-append','-dpsc2',sprintf(fullfile(cwd,'Projects','RS','Networks', ...
%        spm_str_manip(NetFlags.Group,'rt'),[NetFlags.project  '.ps'])));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            %
%   precompute FSLNets       %
%                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath (sprintf('%s/etc/matlab',getenv('FSLDIR')))

addpath (fullfile(cwd,'/matfiles/FSLNets'))
addpath (fullfile(cwd,'/matfiles/FSLNets/L1precision'))
addpath (fullfile(cwd,'/matfiles/FSLNets/pwling'))


ts=nets_load(fullfile(cwd,'Projects','RS','Networks',spm_str_manip(NetFlags.Group,'rt'), ...
    [NetFlags.project '_' des],'TS','TS_all'),tr,1);       % load time series
ts.DD=1:numel(NetFlags.Label);

%   netmat1=nets_netmats(ts,0,'corr');   netmat1=atanh(netmat1);    % full correlation -- with z transformation

netmat2=nets_netmats(ts,0,'icov');   netmat2=atanh(netmat2);    % partial correlation -- with z transformation
netmat3=nets_netmats(ts,0,'ridgep'); netmat3=atanh(netmat3);    % L2-norm Ridge Regression partial correlation
                                                                %   rho=0.1 -with z transformation

%%%     netmats to z
% for i=1:numel(NetFlags.Label)
%     z(:,1:numel(NetFlags.Label),i) = ...
%         netmat1(:,(1+numel(NetFlags.Label)*(i-1)):(numel(NetFlags.Label)+numel(NetFlags.Label)*(i-1)));
% end

for i=1:numel(NetFlags.Label)
    pz(:,1:numel(NetFlags.Label),i) =  ...
        netmat2(:,(1+numel(NetFlags.Label)*(i-1)):(numel(NetFlags.Label)+numel(NetFlags.Label)*(i-1)));
end

for i=1:numel(NetFlags.Label)
    rpz(:,1:numel(NetFlags.Label),i) =  ...
        netmat3(:,(1+numel(NetFlags.Label)*(i-1)):(numel(NetFlags.Label)+numel(NetFlags.Label)*(i-1)));
end

%%%

%%% set outliers (>2*std) to NaN

% for i=1:ana
%     for ii=1:ana
%         z(find(abs(z(:,i,ii))>mean(z(:,i,ii))+2*std(z(:,i,ii))),i,ii)=NaN;
%         zr(find(abs(zr(:,i,ii))>mean(zr(:,i,ii))+2*std(zr(:,i,ii))),i,ii)=NaN;
%         pz(find(abs(pz(:,i,ii))>mean(pz(:,i,ii))+2*std(pz(:,i,ii))),i,ii)=NaN;
%         rpz(find(abs(rpz(:,i,ii))>mean(rpz(:,i,ii))+2*std(rpz(:,i,ii))),i,ii)=NaN;
%     end
% end
% 
% fprintf(1,'%s\n',['For Pearson correlation ' int2str((sum(isnan(z(:)))/numel(z))*100) '% outlieres were excluded (NaN)']);
% fprintf(1,'%s\n',['For Spearman correlation ' int2str((sum(isnan(zr(:)))/numel(zr))*100) '% outlieres were excluded (NaN)']);
% fprintf(1,'%s\n',['For Partial correlation ' int2str((sum(isnan(pz(:)))/numel(z))*100) '% outlieres were excluded (NaN)']);
% fprintf(1,'%s\n',['For regularized Partial correlation ' int2str((sum(isnan(rpz(:)))/numel(zr))*100) '% outlieres were excluded (NaN)']);

NetFlags.z = z;
NetFlags.rawz = z;


NetFlags.origCorr = z;
NetFlags.origRanC = zr;
NetFlags.origParC = pz;
NetFlags.origRegPC = rpz;

Label = NetFlags.Label';
Connectivity = squeeze(mean(z));

save(fullfile(cwd,'Projects','RS','Networks',spm_str_manip(NetFlags.Group,'rt'),[NetFlags.project '_' des], ...
        [NetFlags.project '_' des '.mat']),'z','zr','pz','rpz','Label','Connectivity','tr');

%%% diagnosis

if any(NetFlags.isPat)>0
    groups={'con' 'pat'};
    ts1=nets_load(fullfile(cwd,'Projects','RS','Networks',spm_str_manip(NetFlags.Group,'rt'), ...
        [NetFlags.project '_' des],'TS','TS_con'),tr,1);       % load time series
    ts1.DD=1:numel(NetFlags.Label);
    ts2=nets_load(fullfile(cwd,'Projects','RS','Networks',spm_str_manip(NetFlags.Group,'rt'), ...
        [NetFlags.project '_' des],'TS','TS_pat'),tr,1);       % load time series
    ts2.DD=1:numel(NetFlags.Label);
else
    groups={'all'};
end

%%%% time series diagnostics %%%%

for s=1:numel(groups)
    if numel(groups)==2
        if s==1
            ts=ts1;
        else
            ts=ts2;
        end
    end
    
    %%%%   outlier detection
    [outlier_nodes,outlier_subjects] = nets_outliers(ts);
    if sum(outlier_nodes)~=0 || sum(outlier_subjects)~=0
        fid = fopen(fullfile(cwd,'Projects','RS','Networks',spm_str_manip(NetFlags.Group,'rt'), ...
            [NetFlags.project '_' des],['outlier_' groups{s}  '.txt']),'w+');
        if sum(outlier_nodes)~=0
            outl=find(outlier_nodes);
            for i=1:numel(outl);
                fprintf(fid,'%s\n',['Node #' int2str(outl(i)) ': "' NetFlags.Label{outl(i)} ...
                    '" seems to been an outlier']);
            end
        end
        if sum(outlier_subjects)~=0
            outl=find(outlier_subjects);
            for i=1:numel(outl);
                fprintf(fid,'%s\n',[ 'Subject #' int2str(outl(i)) ': "' NetFlags.sub{outl(i)} ...
                    '" seems to been an outlier']);
            end
        end
        fclose(fid);
    end
    %%%%

    try %%%% produce subject-averaged power spectra 
        [ts_spectra] = nets_spectra_MS(ts);
        movefile(fullfile(cwd,'Projects','RS','Networks',spm_str_manip(NetFlags.Group,'rt'),'FIG.png'),...
            fullfile(cwd,'Projects','RS','Networks',spm_str_manip(NetFlags.Group,'rt'), ...
                [NetFlags.project '_' des],['TS_spectra_' NetFlags.project '_' groups{s} '_' des '.png']));
    end
    
    %%%% TS standard deviation per Node
    netmat0=nets_netmats(ts,0,'amp');               
    [Znet0,Mnet0]=nets_groupmean_amp(netmat0,1);
    movefile(fullfile(cwd,'Projects','RS','Networks',spm_str_manip(NetFlags.Group,'rt'),'FIG.png'),...
        fullfile(cwd,'Projects','RS','Networks',spm_str_manip(NetFlags.Group,'rt'),[NetFlags.project '_' des] ...
            ,['TS_std_' NetFlags.project '_' groups{s} '_' des '.png']));

end

if numel(groups)==2;
    delete(fullfile(cwd,'Projects','RS','Networks',spm_str_manip(NetFlags.Group,'rt'), ...
        [NetFlags.project '_' des],'TS','TS_all','*.*'));
    for s=1:numel(groups)
        copyfile(fullfile(cwd,'Projects','RS','Networks',spm_str_manip(NetFlags.Group,'rt'), ...
            [NetFlags.project '_' des],'TS',['TS_' groups{s}],'*.*'), ...
        fullfile(cwd,'Projects','RS','Networks',spm_str_manip(NetFlags.Group,'rt'), ...
            [NetFlags.project '_' des],'TS','TS_all'))
    end
end

%%%% create  VOI-thumbnails %%%%

if NetFlags.Coordinates==1
   
    NetFlags.X
    NetFlags.Y
    NetFlags.Z
    NetFlags.Radius
end

try
if numel(dir(fullfile(cwd,'VOIs','VOIfiles','thumbnails',[NetFlags.project '.sum'],'*.png'))) ... 
        ~= numel(NetFlags.Label)
    x{numel(NetFlags.Label),1}='';
    if NetFlags.Coordinates==0
        for i=1:numel(NetFlags.VOI);
            x{i,:}=fullfile([NetFlags.VOI{i} ' ']);
        end
    else
        for i=1:numel(NetFlags.Label);
            system([pwd '/createSpheres.sh ' num2str(NetFlags.X(i)) ' ' num2str(NetFlags.Y(i))...
                ' ' num2str(NetFlags.Z(i)) ' ' num2str(NetFlags.Radius) ' ' strrep(NetFlags.Label{i},' ','_')]);
            x{i,:}=fullfile(cwd,'VOIs','VOIfiles','thumbnails',[strrep(NetFlags.Label{i},' ','_') ' ']);
        end
    end
    system([getenv('FSLDIR') '/bin/fslmerge -t ' fullfile(cwd,'VOIs','VOIfiles', ...
        [NetFlags.project '.nii.gz']) ' ' x{:}]);
    try; delete(fullfile(cwd,'VOIs','VOIfiles','thumbnails','*.nii.gz')); end
    system([fullfile(cwd,'matfiles','FSLNets','slices_summary_bin') ' '  fullfile(cwd,'VOIs','VOIfiles', ...
        [NetFlags.project '.nii.gz']) ' 0.1 ' fullfile(cwd,'MaskenEtc','FSLtemplates','MNI152_T1_2mm.nii') ...
            ' ' fullfile(cwd,'VOIs','VOIfiles','thumbnails',[NetFlags.project '.sum']) ' -1']);
    fils = dir(fullfile(cwd,'VOIs','VOIfiles','thumbnails',[NetFlags.project '.sum'])); fils = fils(~[fils.isdir]);
    for file = 1:size(fils,1)
        if numel(strfind(computer,'MAC'))>0
            system(['sips -f horizontal ' fullfile(cwd,'VOIs','VOIfiles','thumbnails', ...
                [NetFlags.project '.sum'],fils(file).name)]);
        else
            system(['convert ' fullfile(cwd,'VOIs','VOIfiles','thumbnails',[NetFlags.project '.sum'], ...
                fils(file).name) ' -flop ' fullfile(cwd,'VOIs','VOIfiles','thumbnails', ...
                    [NetFlags.project '.sum'],fils(file).name)]);     % for Linux
        end
    end
end
end
    
%%%% VOI-thumbnails created %%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            %
%    FSLNets precomputed     %
%                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%  txt output: functional connectivities of all subjects & nodes

wrz = [];
wrzr = [];
wrzp = [];
wrzrp = [];

wn  = {};
cnt = 1;

for i=1:numel(Label)-1
    for ii=i+1:numel(Label)
        wrz(cnt,:) = z(:,i,ii)';
        wrzr(cnt,:) = zr(:,i,ii)';
        wrzp(cnt,:) = pz(:,i,ii)';
        wrzrp(cnt,:) = rpz(:,i,ii)';
        wn{cnt} = [Label{i} ' <-> ' Label{ii}];
        cnt = cnt+1;
    end
end
wrz = wrz';
wrzr = wrzr';
wrzp = wrzp';
wrzrp = wrzrp';

fid = fopen(fullfile(cwd,'Projects','RS','Networks',spm_str_manip(NetFlags.Group,'rt'), ...
    [NetFlags.project '_' des],[NetFlags.project '_' des '_pearson_fc.txt']),'w+');
fprintf(fid,'%s\t','Subject');
fid2 = fopen(fullfile(cwd,'Projects','RS','Networks',spm_str_manip(NetFlags.Group,'rt') ...
    ,[NetFlags.project '_' des],[NetFlags.project '_' des '_spearmann_fc.txt']),'w+');
fprintf(fid2,'%s\t','Subject');
fid3 = fopen(fullfile(cwd,'Projects','RS','Networks',spm_str_manip(NetFlags.Group,'rt'), ...
    [NetFlags.project '_' des],[NetFlags.project '_' des '_partial_fc.txt']),'w+');
fprintf(fid2,'%s\t','Subject');
fid4 = fopen(fullfile(cwd,'Projects','RS','Networks',spm_str_manip(NetFlags.Group,'rt'), ...
    [NetFlags.project '_' des],[NetFlags.project '_' des '_regpar_fc.txt']),'w+');
fprintf(fid2,'%s\t','Subject');


for i=1:numel(wn)
    fprintf(fid,'%s\t',wn{i});
    fprintf(fid2,'%s\t',wn{i});
    fprintf(fid3,'%s\t',wn{i});
    fprintf(fid4,'%s\t',wn{i});
end
for ii=1:numel(NetFlags.sub)
    fprintf(fid,'\n%s\t',NetFlags.sub{ii});
    fprintf(fid2,'\n%s\t',NetFlags.sub{ii});
    fprintf(fid3,'\n%s\t',NetFlags.sub{ii});
    fprintf(fid4,'\n%s\t',NetFlags.sub{ii});
  
    for i=1:numel(wn)
        fprintf(fid,'%8.6f\t',wrz(ii,i));
        fprintf(fid2,'%8.6f\t',wrzr(ii,i));
        fprintf(fid3,'%8.6f\t',wrzp(ii,i));
        fprintf(fid4,'%8.6f\t',wrzrp(ii,i));
    end
end
fclose(fid);
fclose(fid2);
fclose(fid3);
fclose(fid4);
end
