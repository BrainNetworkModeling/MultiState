function fh_FSLNets(des)

warning off

global NetFlags
global fg
global tr
cwd=pwd;

try; NetFlags.gotMask = get(NetFlags.Mask,'value'); end
try; NetFlags.gotPCA = get(NetFlags.PCA,'value'); end
try; NetFlags.gotGlob = get(NetFlags.Glob,'value'); end
try; NetFlags.gotCORR = get(NetFlags.CORR,'value'); end

addpath (sprintf('%s/etc/matlab',getenv('FSLDIR')))
addpath (sprintf('%s/bin',getenv('FSLDIR')))
addpath (sprintf('%s/etc/fslconf',getenv('FSLDIR')))
addpath (fullfile(cwd,'/matfiles/FSLNets'))
addpath (fullfile(cwd,'/matfiles/FSLNets/L1precision'))
addpath (fullfile(cwd,'/matfiles/FSLNets/pwling'))

RSDir = fullfile(cwd,'Projects','RS','Networks');

if nargin ==0
    GM   = {'NoGMmask','IndGMmask','GroupGMmask'};
    PCA  = {'PCA','NoPCA','FIX'};
    Glob = {'TSR','lTSR','WMCSF','lWMCSF','GSR','lGSR','NoGSR'};
    Filt = {'BP','HP'};
    if NetFlags.Coordinates==0
        des = [GM{get(NetFlags.Mask,'value')} '_' PCA{get(NetFlags.PCA,'value')} '_' Glob{get(NetFlags.Glob,'value')} '_' ...
            Filt{get(NetFlags.Filt,'value')}];
    else
        des = [GM{get(NetFlags.Mask,'value')} '_' PCA{get(NetFlags.PCA,'value')}  '_' Glob{get(NetFlags.Glob,'value')} '_' ...
            Filt{get(NetFlags.Filt,'value')} '_' strrep(num2str(NetFlags.Radius,'%2.1f'),'.','')];
    end
end

fprintf(1,'%s\n',['loading time series ...']);
%%% load time-series 
if any(NetFlags.isPat)>0
    groups={'con' 'pat'};
    ts1=nets_load(fullfile(RSDir,spm_str_manip(NetFlags.Group,'rt'),[NetFlags.project '_' des],'TS','TS_con'),tr,1);
    ts1.DD=1:numel(NetFlags.Label);
    ts2=nets_load(fullfile(RSDir,spm_str_manip(NetFlags.Group,'rt'),[NetFlags.project '_' des],'TS','TS_pat'),tr,1);
    ts2.DD=1:numel(NetFlags.Label);
end
ts=nets_load(fullfile(RSDir,spm_str_manip(NetFlags.Group,'rt'),[NetFlags.project '_' des],'TS','TS_all'),tr,1);
ts.DD=1:numel(NetFlags.Label);

ana = {'pearson correlations','spearman rank correlations','partial correlations','regularized parcial correlations (rho=0.1)'};
fprintf(1,'%s\n',['Computing ' ana{NetFlags.gotCORR} ' ...']);

%%% Compute correlations
CORR={'corr';'rank';'icov';'ridgep'};
if any(NetFlags.isPat)>0
    netmat1=nets_netmats(ts1,0,CORR{NetFlags.gotCORR});   netmat1=atanh(netmat1);
    netmat2=nets_netmats(ts2,0,CORR{NetFlags.gotCORR});   netmat2=atanh(netmat2);
end
netmat=nets_netmats(ts,0,CORR{NetFlags.gotCORR});   netmat=atanh(netmat);

thumbDir = fullfile(cwd,'VOIs','VOIfiles','thumbnails');
%%%% create  VOI-thumbnails %%%%
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
            system([cwd '/createSpheres.sh ' num2str(NetFlags.X(i)) ' ' num2str(NetFlags.Y(i))...
                ' ' num2str(NetFlags.Z(i)) ' ' num2str(NetFlags.Radius) ' ' strrep(NetFlags.Label{i},' ','_')]);
            x{i,:}=fullfile(cwd,'VOIs',[strrep(NetFlags.Label{i},' ','_') ' ']);
        end
    end
    system([getenv('FSLDIR') '/bin/fslmerge -t ' fullfile(cwd,'VOIs','VOIfiles', ...
        [NetFlags.project '.nii.gz']) ' ' x{:}]);
%     try; delete(fullfile(cwd,'VOIs','VOIfiles','thumbnails','*.nii.gz')); end
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
%     if numel(dir(fullfile(thumbDir,[NetFlags.project '.sum'],'*.png'))) ~= numel(NetFlags.VOI)
%         fprintf(1,'%s\n','Creating VOI thumbnail slices ...');
%         x{numel(NetFlags.VOI),1}='';
%         for i=1:numel(NetFlags.VOI); x{i,:}=fullfile([NetFlags.VOI{i} ' ']); end
%         system([getenv('FSLDIR') '/bin/fslmerge -t ' fullfile(cwd,'VOIs','VOIfiles',[NetFlags.project '.nii.gz']) ' ' x{:}]);
%         system([fullfile(cwd,'matfiles','FSLNets','slices_summary_bin') ' ' fullfile(cwd,'VOIs','VOIfiles', ...
%             [NetFlags.project '.nii.gz']) ' 0.1 ' fullfile(cwd,'MaskenEtc','FSLtemplates','MNI152_T1_2mm.nii') ' ' ...
%             fullfile(thumbDir,[NetFlags.project '.sum']) ' -1']);
%         fils = dir(fullfile(thumbDir,[NetFlags.project '.sum'])); fils = fils(~[fils.isdir]);
%         for file = 1:size(fils,1)
%             if numel(strfind(computer,'MAC'))>0
%                 system(['sips -f horizontal ' fullfile(thumbDir,[NetFlags.project '.sum'], fils(file).name)]);
%             else
%                 system(['convert ' fullfile(thumbDir,[NetFlags.project '.sum'],fils(file).name) ...
%                     ' -flop ' fullfile(thumbDir,[NetFlags.project '.sum'],fils(file).name)]);
%             end
%         end
%     end
end

%%% print analysis settings
print(fg,'-append','-dpsc2',sprintf(fullfile(RSDir,spm_str_manip(NetFlags.Group,'rt'), ...
    [NetFlags.project '_' des],'clustering.ps')));

fprintf(1,'%s\n','Computing group mean connectivity matrices & hierarchical clustering ...');
%%% estimate group mean/one-group-t-test and consistency of netmats across runs/subjects
if any(NetFlags.isPat)>0
    [Znet1,Mnet1]=nets_groupmean_MS(netmat1,1,des);
    [Znet2,Mnet2]=nets_groupmean_MS(netmat2,1,des);
    if sum(NetFlags.isPat==1)==sum(NetFlags.isPat==0)
        [Znet,Mnet]=nets_groupmean_MS(netmat,0,des,sum(NetFlags.isPat==1));
    else
        [Znet,Mnet]=nets_groupmean_MS(netmat,0,des);
    end
else
    [Znet,Mnet]=nets_groupmean_MS(netmat,1,des);
end


%%% hierarchical clustering of crosscorrelation
if any(NetFlags.isPat)>0
    nets_hierarchy_MS(Znet1,Mnet1,1:numel(NetFlags.Label),fullfile(thumbDir,[NetFlags.project '.sum']),des);
    nets_hierarchy_MS(Znet2,Mnet2,1:numel(NetFlags.Label),fullfile(thumbDir,[NetFlags.project '.sum']),des);
    nets_netweb(Znet1,Mnet1,1:numel(NetFlags.Label),fullfile(thumbDir,[NetFlags.project]),...
        fullfile(RSDir,spm_str_manip(NetFlags.Group,'rt'),[NetFlags.project '_' des],'netweb_Con'))
    nets_netweb(Znet2,Mnet2,1:numel(NetFlags.Label),fullfile(thumbDir,[NetFlags.project]),...
        fullfile(RSDir,spm_str_manip(NetFlags.Group,'rt'),[NetFlags.project '_' des],'netweb_Pat'))
else
    nets_hierarchy_MS(Znet,Mnet,1:numel(NetFlags.Label),fullfile(thumbDir,[NetFlags.project '.sum']),des);
    nets_netweb(Znet,Mnet,1:numel(NetFlags.Label),fullfile(thumbDir,[NetFlags.project]),...
        fullfile(RSDir,spm_str_manip(NetFlags.Group,'rt'),[NetFlags.project '_' des],'netweb'))
end



%     if numel(NetFlags.Label)>=6
%         showN = 15;
%     elseif numel(NetFlags.Label)==5
%         showN = 10;
%     elseif numel(NetFlags.Label)==4
%         showN = 6;
%     elseif numel(NetFlags.Label)==3
%         showN = 3;
%     end
%     nets_edgepics_MS(1:numel(NetFlags.Label),fullfile(thumbDir,[NetFlags.project '.sum']),Mnet,Znet,showN,des);


% if any(NetFlags.isPat)>0
%     [p_uncorrected1,p_corrected1] = nets_glm_MS(netmat,fullfile(cwd,'matfiles','FSLNets',[num2str(sum(NetFlags.isPat==0)) ...
%         '_design.mat']),fullfile(cwd,'matfiles','FSLNets',[num2str(sum(NetFlags.isPat==0)) '_design.con']),1,des);
% end

end