function fh_StrucCov(varargin)

%%% structural covariance %%%

global NetFlags
cwd=pwd;

addpath (sprintf('%s/etc/matlab',getenv('FSLDIR')))

addpath (fullfile(cwd,'/matfiles/FSLNets'))
addpath (fullfile(cwd,'/matfiles/FSLNets/L1precision'))
addpath (fullfile(cwd,'/matfiles/FSLNets/pwling'))

VBMDir = fullfile(cwd,'Projects','VBM','Regional');
thumbDir = fullfile(cwd,'VOIs','VOIfiles','thumbnails');
%%%% create  VOI-thumbnails %%%%

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



switch get(NetFlags.Corr,'value')
    case 1
        [r,p] = corr(NetFlags.z(:,:));
    case 2
        [r,p] = corr(NetFlags.z(:,:),'type','Spearman');
    case 3
        [r,p] = partialcorr(NetFlags.z);
end

r(r==1)=0;
NetFlags.c = r;

%%% Threshold %%%

toDo = numel(NetFlags.Label);
toDo = (toDo*toDo-toDo)/toDo;

switch get(NetFlags.Method,'value')
    case 1
        pThr    = NetFlags.pvalue;
        desc    = ['p < ' num2str(NetFlags.pvalue) ' (uncorr.)'];
    case 2
        pThr    = NetFlags.pvalue/toDo;
        desc    = ['p < ' num2str(NetFlags.pvalue) ' (Bonf. corr.)'];
    case 3
        pp = triu(p,1); pp = pp(:); pp(pp==0) = [];
        pThr    = max(spm_uc_FDR(NetFlags.pvalue,1,'P',1,sort(pp,'ascend'),[]),NetFlags.pvalue/toDo);
        desc    = ['p < ' num2str(NetFlags.pvalue) ' (FDR corr.)'];
end


p(p==0)=1; p=1-p(:,:); p(p<=1-pThr)=0;


Action = varargin{1};
switch lower(Action),
    case 'corr'

        % correlations %
        pc=p; pc(find(pc(:)>0 & r(:)<0))=0;
        showN = sum(pc(:)>0)/2; if showN >= 35; showN = 35; end
        if showN > 0
            nets_edgepics_MS2(1:numel(NetFlags.Label),fullfile(thumbDir,[NetFlags.project '.sum']),r,pc,showN);
        end

        % anti-correlations %
        pac=p; pac(find(pac(:)>0 & r(:)>0))=0;
        showN = sum(pac(:)>0)/2; if showN >= 35; showN = 35; end
        if showN > 0
            nets_edgepics_MS2(1:numel(NetFlags.Label),fullfile(thumbDir,[NetFlags.project '.sum']),r,pac,showN);
        end
   
    case 'ward'
        % Ward clustering of structural covariance
        nets_hierarchy_MS2(r,p,1:numel(NetFlags.Label),fullfile(thumbDir,[NetFlags.project '.sum']));
        nets_netweb(r,p,1:numel(NetFlags.Label),fullfile(thumbDir,[NetFlags.project]),...
            fullfile(VBMDir,spm_str_manip(NetFlags.Group,'rt'),'netweb'))

end


% fprintf('\n%s \t', '&:-{D');




% z=[];
% z(:,:) = atanh(r);
% zr(:,:) = atanh(rr);  
% z(z==Inf)=0;
% zr(zr==Inf)=0;

% [Znet,Mnet]=nets_groupmean(z,1);
    


