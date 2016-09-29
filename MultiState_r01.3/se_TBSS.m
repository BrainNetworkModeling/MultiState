function se_TBSS(varargin)


global GLMFlags
clc


if (nargin==0), Action = 'init'; else, Action = varargin{1}; end

switch lower(Action),
    
    case 'init'
        fg = spm_figure('GetWin','Graphics');
        Finter = spm('CreateIntWin','on');
        spm_figure('Clear','Graphics');
        FS   = spm('FontSizes');
        
        ok = 1;
        if ~isfield(GLMFlags,'sub')
            GLMFlags.Group = 'No group lookup selected';
            ok = 0;
        end
        uicontrol(fg,'Style','Pushbutton','Units','normalized','Position',[0.05 .95 .5 .04],'String',{spm_str_manip(GLMFlags.Group,'rt')},...
            'HorizontalAlignment','center','FontWeight','bold','Callback','se_TBSS(''gselect'')','FontSize',10);
        
        if ~isfield(GLMFlags,'Label')
            GLMFlags.project = 'No VOI file selected';
            GLMFlags.Coordinates = 0;
        end
        
        uicontrol(fg,'Style','Pushbutton','Units','normalized','Position',[0.1 .9 .4 .04],'String',{spm_str_manip(GLMFlags.project)},...
            'HorizontalAlignment','center','FontWeight','bold','Callback','se_TBSS(''vselect'')','FontSize',10);

        
        
        if ok == 0
            uicontrol(fg,'style','pushbutton','units','normalized','position',[0.65 .9 .3 .09],'string',{'Compute'},...
                'horizontalalignment','center','fontweight','bold','foregroundcolor','w','FontSize',16);
        else
            uicontrol(fg,'style','pushbutton','units','normalized','position',[0.65 .9 .3 .09],'string',{'Compute'},...
                'horizontalalignment','center','fontweight','bold','callback','se_computeTBSS','foregroundcolor','r','FontSize',16);
        end
        
        
        
        
        
        if ok == 1

                try; sel = get(GLMFlags.Analysis,'Value'); catch; sel = 1; end
                GLMFlags.Analysis = uicontrol(fg,'Style','popupmenu','Units','normalized','Position',[0.02 .65 .63 .1],...
                    'String',str2mat('Group differences','Group differences with covariate','Group x Covariate',...
                    'All: Covariate Analysis','Patients: Covariate Analysis'),'Value',sel,'ToolTipString','Select analysis model',...
                    'FontSize',14,'Callback','se_TBSS(''amode'')');
                GLMFlags.HasCov = [0 1 1 1 1];
                
                
                
            try; sel = get(GLMFlags.Conf1,'Value'); catch; sel = 1; end
            GLMFlags.Conf1 = uicontrol(fg,'Style','popupmenu','Units','normalized','Position',[0.01 .76 .24 .07],...
                'String',{'Confound 1: None' GLMFlags.Covariates{:}},'Value',sel,'ToolTipString','Select confound factor',...
                'FontSize',8);
            try; sel = get(GLMFlags.Conf2,'Value'); catch; sel = 1; end
            GLMFlags.Conf2 = uicontrol(fg,'Style','popupmenu','Units','normalized','Position',[0.25 .76 .24 .07],...
                'String',{'Confound 2: None' GLMFlags.Covariates{:}},'Value',sel,'ToolTipString','Select confound factor',...
                'FontSize',8);
            try; sel = get(GLMFlags.Conf3,'Value'); catch; sel = 1; end
            GLMFlags.Conf3 = uicontrol(fg,'Style','popupmenu','Units','normalized','Position',[0.50 .76 .24 .07],...
                'String',{'Confound 3: None' GLMFlags.Covariates{:}},'Value',sel,'ToolTipString','Select confound factor',...
                'FontSize',8);
            try; sel = get(GLMFlags.Confx,'Value'); catch; sel = 1; end
            GLMFlags.Confx = uicontrol(fg,'Style','popupmenu','Units','normalized','Position',[0.75 .76 .24 .07],...
                'String',{'No removal','1st Order','2nd Order','3rd Order'},'Value',sel,'ToolTipString','Select to remove confounds',...
                'FontSize',8,'foregroundcolor','r');
                
            if isfield(GLMFlags,'Label')
                try; sel = get(NetFlags.VBM,'Value'); catch; sel = 1; end
                GLMFlags.VBM = uicontrol(fg,'Style','popupmenu','Units','normalized','Position',[0.15 .78 .3 .1],...
                    'String',str2mat('','m Volume as Covariate','m0 Volume as Covariate'),'Value',sel,'FontSize',8,'ToolTipString','Remove global signal drifts');
            end
                
        end
        
        
        uicontrol(fg,'Style','PushButton','Units','normalized','Position',[.9 .02 .08 .05],'ForegroundColor','r','FontWeight','bold',...
            'String','EXIT','Callback','se_TBSS(''exit'');','ToolTipString','quit');
        
    case 'amode'
        fg = spm_figure('GetWin','Graphics');
        if GLMFlags.HasCov(get(GLMFlags.Analysis,'Value'))>0
            GLMFlags.ActiveCov = [];
            GLMFlags.Cselect = uicontrol(fg,'Style','popupmenu','Units','normalized','Position',[0.02 .6 .63 .1],...
                'String',GLMFlags.Covariates,'Value',1,'ToolTipString','Select covariates',...
                'FontSize',14,'Callback','se_TBSS(''csel'')');
            GLMFlags.Cdel1 = uicontrol(fg,'Style','PushButton','Units','normalized','Position',[0.68 .66 .3 .04],...
                'String','Delete last covariate','FontSize',14,'Callback','se_TBSS(''cdel1'')');
            GLMFlags.Cdel0 = uicontrol(fg,'Style','PushButton','Units','normalized','Position',[0.68 .71 .3 .04],...
                'String','Delete  ALL  covariates','FontSize',14,'Callback','se_TBSS(''cdel0'')');

%             GLMFlags.ActiveNui = [];
%             GLMFlags.Nselect = uicontrol(fg,'Style','popupmenu','Units','normalized','Position',[0.02 .3 .63 .1],...
%                 'String',GLMFlags.Covariates,'Value',1,'ToolTipString','Select confounds',...
%                 'FontSize',14,'Callback','se_TBSS(''nsel'')');
%             GLMFlags.Ndel1 = uicontrol(fg,'Style','PushButton','Units','normalized','Position',[0.68 .36 .3 .04],...
%                 'String','Delete last confound','FontSize',14,'Callback','se_TBSS(''ndel1'')');
%             GLMFlags.Ndel0 = uicontrol(fg,'Style','PushButton','Units','normalized','Position',[0.68 .31 .3 .04],...
%                 'String','Delete  ALL  confounds','FontSize',14,'Callback','se_TBSS(''ndel0'')');

            
        else
            try; delete(GLMFlags.Cselect); end
            try; delete(GLMFlags.Cdel1); end
            try; delete(GLMFlags.Cdel0); end
            try; delete(GLMFlags.CovNames); end
            
        end
        
    case 'csel'
        fg = spm_figure('GetWin','Graphics');
        if ~any(GLMFlags.ActiveCov == get(GLMFlags.Cselect,'Value'))
            GLMFlags.ActiveCov = [GLMFlags.ActiveCov get(GLMFlags.Cselect,'Value')];
            GLMFlags.CovNames  = uicontrol(fg,'Style','text','Units','normalized','Position',[0.02 .4 .63 .2],...
                'HorizontalAlignment','center','BackgroundColor',[1 1 1],...
                'String',{GLMFlags.Covariates{GLMFlags.ActiveCov}},'FontSize',14);
        end
        
    case 'cdel1'
        fg = spm_figure('GetWin','Graphics');
        if numel(GLMFlags.ActiveCov>0)
            GLMFlags.ActiveCov = GLMFlags.ActiveCov(1:end-1);
            GLMFlags.CovNames  = uicontrol(fg,'Style','text','Units','normalized','Position',[0.02 .4 .63 .2],...
                'HorizontalAlignment','center','BackgroundColor',[1 1 1],...
                'String',{GLMFlags.Covariates{GLMFlags.ActiveCov}},'FontSize',14);
        end
        
    case 'cdel0'
        fg = spm_figure('GetWin','Graphics');
        if numel(GLMFlags.ActiveCov>0)
            GLMFlags.ActiveCov = [];
            GLMFlags.CovNames  = uicontrol(fg,'Style','text','Units','normalized','Position',[0.02 .4 .63 .2],...
                'HorizontalAlignment','center','BackgroundColor',[1 1 1],...
                'String',{GLMFlags.Covariates{GLMFlags.ActiveCov}},'FontSize',14);
        end
        
        
        
        
        
        
        
        
    case 'gselect'
        GLMFlags.Group = spm_select(1,'any',['Select group description'],[],fullfile(pwd,'Lookups'),'.mat',1);
        tmp = load(GLMFlags.Group);
        GLMFlags.Covariates = tmp.Covariates;
        GLMFlags.sub = tmp.sub(tmp.hasDTI==1);
        GLMFlags.Cov = tmp.Cov(tmp.hasDTI==1,:);
        GLMFlags.rawCov = tmp.Cov(tmp.hasDTI==1,:);
        GLMFlags.isPat = tmp.isPat(tmp.hasDTI==1);
        GLMFlags.SubDir = tmp.SubDir(tmp.hasDTI==1);
        GLMFlags.Covariates = tmp.Covariates;
        GLMFlags.rawCovariates = tmp.Covariates;
        if size(GLMFlags.isPat,1)>size(GLMFlags.isPat,2)
            GLMFlags.isPat = GLMFlags.isPat';
        end
        clear tmp
        
%         if numel(GLMFlags.Covariates)>0
%             a = get(0,'ScreenSize');
%             FS   = spm('FontSizes');
%             PF   = spm_platform('fonts');
%             try
%                 F = findall(allchild(0),'Flat','Tag','CovBrowser');
%                 GLMFlags = rmfield(GLMFlags,'Cfg');
%                 delete(F)
%             end
%             GLMFlags.Cfg =  figure(...
%                 'Tag','CovBrowser',                          'Position',[a(3)*.01 a(4)*.3 a(3)*.4 a(4)*.6],...
%                 'Resize','on',                            'MenuBar','figure',                 'Color','w',                              'ColorMap',gray(64),...
%                 'DefaultTextColor','k',                   'DefaultTextInterpreter','tex',     'DefaultTextFontName',PF.helvetica,       'DefaultTextFontSize',FS(12),...
%                 'DefaultAxesColor','w',                   'DefaultAxesXColor','k',            'DefaultAxesYColor','k',                  'DefaultAxesZColor','k',...
%                 'DefaultAxesFontName',PF.helvetica,       'DefaultPatchFaceColor','k',        'DefaultPatchEdgeColor','k',              'DefaultSurfaceEdgeColor','k',...
%                 'DefaultLineColor','k',                   'DefaultUicontrolFontName',PF.helvetica,...
%                 'DefaultUicontrolFontSize',FS(12),        'DefaultUicontrolInterruptible','on', 'PaperType','A4',                         'PaperUnits','normalized',...
%                 'PaperPosition',[.0726 .0644 .854 .870],  'InvertHardcopy','off',             'Renderer','zbuffer',                     'Visible','on','Name','Results');
%             GLMFlags.Cfg.select = uicontrol(GLMFlags.Cfg,'Style','popupmenu','Units','normalized','Position',[0.05 .9 .9 .1],...
%                 'String',GLMFlags.Covariates,'Value',1,'ToolTipString','Select covariates',...
%                 'FontSize',14,'Callback','se_TBSS(''cdisp'')');
%         else
%             try
%                 F = findall(allchild(0),'Flat','Tag','CovBrowser');
%                 GLMFlags = rmfield(GLMFlags,'Cfg');
%                 delete(F)
%             end
%         end
        se_TBSS('init')
        
        
    case 'cdisp'
        F = findall(allchild(0),'Flat','Tag','CovBrowser');
        selected = get(GLMFlags.Cfg.select,'Value');
        try
            delete(GLMFlags.Cfg.show)
        end
        GLMFlags.Cfg.show = axes('Parent',F,'Units','normalized','Position',[0.05 .05 .9 .8]);
        cla
        
        Q = find(~isnan(GLMFlags.Cov(:,selected)));
        
        if any(GLMFlags.isPat)
            fac = 10;
        else
            fac = 5;
        end
        
        if numel(unique(GLMFlags.Cov(Q,selected)))<round(numel((GLMFlags.sub))/fac)
            [n xout] = hist(GLMFlags.Cov(Q,selected),unique(GLMFlags.Cov(Q,selected)));
        else
            if round(numel((GLMFlags.sub))/fac)<5
                [n xout] = hist(GLMFlags.Cov(Q,selected),5);
            elseif round(numel((GLMFlags.sub))/fac)>20
                [n xout] = hist(GLMFlags.Cov(Q,selected),20);
            else
                [n xout] = hist(GLMFlags.Cov(Q,selected),round(numel((GLMFlags.sub))/fac));
            end
        end
        
        
        if any(GLMFlags.isPat)
            [n1 xx] = hist(GLMFlags.Cov(Q(GLMFlags.isPat(Q)==0),selected),xout);
            [n2 xx] = hist(GLMFlags.Cov(Q(GLMFlags.isPat(Q)==1),selected),xout);
            try
                n = [n1; n2]';
            catch
                n = n2;
            end
        end
        
        bar(GLMFlags.Cfg.show,n); axis tight
        if numel(unique(round(xout))) == numel(unique((xout)))
            xout = round(xout);
        else
            xout = round(xout*10)/10;
        end
        set(gca,'XTick',[1:numel(xout)],'XTickLabel',xout)
        
        if any(GLMFlags.isPat)
            tit = {['Control mean: ' num2str(mean(GLMFlags.Cov(Q(GLMFlags.isPat(Q)==0),selected)),'%4.1f')];...
                ['Patient mean: ' num2str(mean(GLMFlags.Cov(Q(GLMFlags.isPat(Q)==1),selected)),'%4.1f')]};
        else
            tit = {['Group mean: ' num2str(mean(GLMFlags.Cov(Q,selected)),'%4.1f')];...
                ['Distinct values: ' int2str(numel(unique(GLMFlags.Cov(Q,selected))))]};
        end
        title(tit)
        
    case 'vselect'
        TxtName  = spm_select(1,'any',['Select coordinate or VOI list'],[],fullfile(pwd,'VOIs'),'txt',1);
        try
            [GLMFlags.X GLMFlags.Y GLMFlags.Z GLMFlags.Label]     = textread(TxtName,'%f %f %f %s');
            GLMFlags.Coordinates = 1;
        catch
            [GLMFlags.VOI GLMFlags.Label] = textread(TxtName,'%s %s');
            GLMFlags.Coordinates = 0;
        end
        GLMFlags.Covariates = {GLMFlags.rawCovariates{:} GLMFlags.Label{:}};
        
        
        GLMFlags.project = strrep(strrep(spm_str_manip(TxtName,'rt'),'_VOIs',''),'_ImageList','');
        se_TBSS('init')
        
        
        
    case 'exit'
        delete(gcf)
%        se_figure('Clear','Graphics');
        clear
        OverGUI;
end


