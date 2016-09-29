function se_GLM(varargin)


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
            'HorizontalAlignment','center','FontWeight','bold','Callback','se_GLM(''gselect'')',...
            'ToolTipString','Choose Lookup-file describing resting-state sample');
        
        if ~isfield(GLMFlags,'Label')
            GLMFlags.project = 'No VOI file selected';
            GLMFlags.Coordinates = 0;
            ok = 0;
        end
        uicontrol(fg,'Style','Pushbutton','Units','normalized','Position',[0.05 .9 .5 .04],'String',{spm_str_manip(GLMFlags.project)},...
            'HorizontalAlignment','center','FontWeight','bold','Callback','se_GLM(''vselect'')',...
            'ToolTipString','File may contain coordinates or image VOIs');
        
        if ok == 0
            uicontrol(fg,'style','pushbutton','units','normalized','position',[0.65 .9 .3 .09],'string',{'Compute'},...
                'horizontalalignment','center','fontweight','bold','foregroundcolor','w');
        else
            uicontrol(fg,'style','pushbutton','units','normalized','position',[0.65 .9 .3 .09],'string',{'Compute'},...
                'horizontalalignment','center','fontweight','bold','callback','se_computeGLM','foregroundcolor','r');

            if numel(strfind(computer,'MAC'))>0
                uicontrol(fg,'Style','PushButton','Units','normalized','Position',[.02 .02 .15 .05],'ForegroundColor','b','FontWeight','bold',...
                    'String','sent GLM','callback','se_remGLM;','ToolTipString','sent GLM & VOIs to FZJ');
            end
 
        end
        
        
        try; sel = get(GLMFlags.Mask,'Value'); catch; sel = 1; end
        GLMFlags.Mask = uicontrol(fg,'Style','popupmenu','Units','normalized','Position',[0.02 .8 .18 .08],...
            'String',str2mat('VOIs: No GM mask','Individual median-split','Group median-split'),'FontSize',8,'Value',sel,...
            'ToolTipString','Select method for VOI gray matter masking');
        
        try; if mean(GLMFlags.fix)>.7
                GLMFlags.PCA = uicontrol(fg,'Style','popupmenu','Units','normalized','Position',[0.21 .8 .18 .08],...
                    'String',str2mat('PCA denoising','No PCA denoising','FIX de-noising'),'Value',sel,'foregroundcolor','r','FontSize',10,...
                    'ToolTipString','FMRIB''s ICA-based X-noiseifier');
            else
                GLMFlags.PCA = uicontrol(fg,'Style','popupmenu','Units','normalized','Position',[0.21 .8 .18 .08],...
                    'String',str2mat('PCA denoising','No PCA denoising'),'Value',sel,'FontSize',8,'ToolTipString','Include 5 component PCA denoising');
            end
        catch
            try; sel = GLMFlags.PCA; catch; sel = 2; end
            GLMFlags.PCA = uicontrol(fg,'Style','popupmenu','Units','normalized','Position',[0.21 .8 .18 .08],...
                'String',str2mat('PCA denoising','No PCA denoising'),'Value',sel,'FontSize',8,'ToolTipString','Include 5 component PCA denoising');
        end 
        
        try; sel = get(GLMFlags.Glob,'Value'); catch; sel = 7; end
        GLMFlags.Glob = uicontrol(fg,'Style','popupmenu','Units','normalized','Position',[0.4 .8 .22 .08],'String',str2mat...
            ('Full tissue signal regression','Linear tissue signal regression','Full WM/CSF signal regression','Linear WM/CSF signal regression',...
            'Full global signal regression','Linear global signal regression', 'No Global signal regression'),...
            'Value',sel,'FontSize',8,'ToolTipString','Remove global signal drifts');
        
        try; sel = get(GLMFlags.Filt,'Value'); catch; sel = 1; end
        GLMFlags.Filt = uicontrol(fg,'Style','popupmenu','Units','normalized','Position',[0.63 .8 .13 .08],...
            'String',str2mat('Band-pass filter','High-pass filter'),'Value',sel,'FontSize',8,'ToolTipString','Filter for BOLD frequencies');
  
        
        
        if GLMFlags.Coordinates>0
            if ~isfield(GLMFlags,'Radius'); GLMFlags.Radius = 5; end
            uicontrol(fg,'Style','Pushbutton','Units','normalized','Position',[0.78 .84 .2 .04],'String',{['Sphere: ' num2str(GLMFlags.Radius) 'mm radius']},...
                'HorizontalAlignment','center','FontWeight','bold','Callback','se_GLM(''radius'')','ForegroundColor','k');
        end
        
        
        if ok == 1
        % Select blocking variable 
        
            GLMFlags.Block = uicontrol(fg,'Style','popupmenu','Units','normalized','Position',[0.02 .76 .35 .07],...
                'String',{'No group blocking factor' GLMFlags.Covariates{:}},'Value',1,'ToolTipString','Select covariates',...
                'FontSize',8);

            
        try; sel = get(GLMFlags.VBM,'Value'); catch; sel = 1; end
        GLMFlags.VBM = uicontrol(fg,'Style','popupmenu','Units','normalized','Position',[0.37 .76 .3 .07],...
            'String',str2mat('','m Volume as Confound','m0 Volume as Confound'),'Value',sel,'FontSize',8,'ToolTipString','Account for atrophy');
            
        
        
        
            try; sel = get(GLMFlags.Conf1,'Value'); catch; sel = 1; end
            GLMFlags.Conf1 = uicontrol(fg,'Style','popupmenu','Units','normalized','Position',[0.01 .71 .24 .07],...
                'String',{'Confound 1: None' GLMFlags.Covariates{:}},'Value',sel,'ToolTipString','Select confound factor',...
                'FontSize',8);
            try; sel = get(GLMFlags.Conf2,'Value'); catch; sel = 1; end
            GLMFlags.Conf2 = uicontrol(fg,'Style','popupmenu','Units','normalized','Position',[0.25 .71 .24 .07],...
                'String',{'Confound 2: None' GLMFlags.Covariates{:}},'Value',sel,'ToolTipString','Select confound factor',...
                'FontSize',8);
            try; sel = get(GLMFlags.Conf3,'Value'); catch; sel = 1; end
            GLMFlags.Conf3 = uicontrol(fg,'Style','popupmenu','Units','normalized','Position',[0.50 .71 .24 .07],...
                'String',{'Confound 3: None' GLMFlags.Covariates{:}},'Value',sel,'ToolTipString','Select confound factor',...
                'FontSize',8);
            try; sel = get(GLMFlags.Confx,'Value'); catch; sel = 1; end
            GLMFlags.Confx = uicontrol(fg,'Style','popupmenu','Units','normalized','Position',[0.75 .71 .24 .07],...
                'String',{'No removal','1st Order','2nd Order','3rd Order'},'Value',sel,'ToolTipString','Select to remove confounds',...
                'FontSize',8,'foregroundcolor','r');
            
            try; sel = get(GLMFlags.Conf1t,'Value'); catch; sel = 1; end
            GLMFlags.Conf1t = uicontrol(fg,'Style','popupmenu','Units','normalized','Position',[0.01 .68 .24 .07],...
                'String',str2mat('Categorical','Continuous'),'Value',sel,'ToolTipString','Select confound factor',...
                'FontSize',8);
            try; sel = get(GLMFlags.Conf2t,'Value'); catch; sel = 1; end
            GLMFlags.Conf2t = uicontrol(fg,'Style','popupmenu','Units','normalized','Position',[0.25 .68 .24 .07],...
                'String',str2mat('Categorical','Continuous'),'Value',sel,'ToolTipString','Select confound factor',...
                'FontSize',8);
            try; sel = get(GLMFlags.Conf3t,'Value'); catch; sel = 1; end
            GLMFlags.Conf3t = uicontrol(fg,'Style','popupmenu','Units','normalized','Position',[0.50 .68 .24 .07],...
                'String',str2mat('Categorical','Continuous'),'Value',sel,'ToolTipString','Select confound factor',...
                'FontSize',8);
            
            
            
            
            if any(GLMFlags.isPat)
                try; sel = get(GLMFlags.Analysis,'Value'); catch; sel = 1; end
                GLMFlags.Analysis = uicontrol(fg,'Style','popupmenu','Units','normalized','Position',[0.02 .6 .63 .1],...
                    'String',str2mat('Group differences','Group differences with covariate','Group x Covariate',...
                    'Patients: Covariates individually','Patients: Multiple Covariates','Paired_t-test','Paired t-test with covariate'),'Value',sel,...
                    'ToolTipString','Select analysis model','FontSize',14,'Callback','se_GLM(''amode'')');
                GLMFlags.HasCov = [0 1 1 1 1 0 1];
            else
                try; sel = get(GLMFlags.Analysis,'Value'); catch; sel = 1; end
                GLMFlags.Analysis = uicontrol(fg,'Style','popupmenu','Units','normalized','Position',[0.02 .6 .63 .1],...
                    'String',str2mat('Main effect','Covariates individually','Multiple Covariates'),'Value',sel,'ToolTipString','Select analysis model',...
                    'FontSize',14,'Callback','se_GLM(''amode'')');
                GLMFlags.HasCov = [0 1 1];
            end
        end
        
        uicontrol(fg,'Style','PushButton','Units','normalized','Position',[.9 .02 .08 .05],'ForegroundColor','r','FontWeight','bold',...
            'String','EXIT','Callback','se_GLM(''exit'');','ToolTipString','quit');
        
    case 'amode'
        fg = spm_figure('GetWin','Graphics');
        if GLMFlags.HasCov(get(GLMFlags.Analysis,'Value'))>0
            GLMFlags.ActiveCov = [];
            GLMFlags.Cselect = uicontrol(fg,'Style','popupmenu','Units','normalized','Position',[0.02 .6 .63 .1],...
                'String',GLMFlags.Covariates,'Value',1,'ToolTipString','Select covariates',...
                'FontSize',14,'Callback','se_GLM(''csel'')');
            GLMFlags.Cdel1 = uicontrol(fg,'Style','PushButton','Units','normalized','Position',[0.68 .6 .3 .04],...
                'String','Delete last covariate','FontSize',14,'Callback','se_GLM(''cdel1'')');
            GLMFlags.Cdel0 = uicontrol(fg,'Style','PushButton','Units','normalized','Position',[0.68 .66 .3 .04],...
                'String','Delete  ALL  covariates','FontSize',14,'Callback','se_GLM(''cdel0'')');
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
            GLMFlags.CovNames  = uicontrol(fg,'Style','text','Units','normalized','Position',[0.02 .2 .53 .4],...
                'HorizontalAlignment','center','BackgroundColor',[1 1 1],...
                'String',{GLMFlags.Covariates{GLMFlags.ActiveCov}},'FontSize',14);
        end
        
    case 'cdel1'
        fg = spm_figure('GetWin','Graphics');
        if numel(GLMFlags.ActiveCov>0)
            GLMFlags.ActiveCov = GLMFlags.ActiveCov(1:end-1);
            GLMFlags.CovNames  = uicontrol(fg,'Style','text','Units','normalized','Position',[0.02 .2 .53 .4],...
                'HorizontalAlignment','center','BackgroundColor',[1 1 1],...
                'String',{GLMFlags.Covariates{GLMFlags.ActiveCov}},'FontSize',14);
        end
        
    case 'cdel0'
        fg = spm_figure('GetWin','Graphics');
        if numel(GLMFlags.ActiveCov>0)
            GLMFlags.ActiveCov = [];
            GLMFlags.CovNames  = uicontrol(fg,'Style','text','Units','normalized','Position',[0.02 .2 .53 .4],...
                'HorizontalAlignment','center','BackgroundColor',[1 1 1],...
                'String',{GLMFlags.Covariates{GLMFlags.ActiveCov}},'FontSize',14);
        end
        
        
    case 'gselect'
        GLMFlags.Group = spm_select(1,'mat',['Select group description'],[],fullfile(pwd,'Lookups'),'.mat',1);
        tmp = load(GLMFlags.Group);
        GLMFlags.Covariates = tmp.Covariates;
        GLMFlags.sub = tmp.sub(tmp.hasRS==1);
        GLMFlags.Cov = tmp.Cov(tmp.hasRS==1,:);
        GLMFlags.isPat = tmp.isPat(tmp.hasRS==1);
        GLMFlags.SubDir = tmp.SubDir(tmp.hasRS==1);
        GLMFlags.Covariates = tmp.Covariates;
        GLMFlags.matchQ=tmp.matchQ;
        GLMFlags.fix = tmp.hasRSfix(tmp.hasRS==1);
        
        if size(GLMFlags.isPat,1)>size(GLMFlags.isPat,2)
            GLMFlags.isPat = GLMFlags.isPat';
        end
        clear tmp
        
        se_GLM('init')
        
        
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
        TxtName  = spm_select(1,'any',['Select coordinate or VOI list'],[],fullfile(pwd,'VOIs'),'.txt',1);
        try
            [GLMFlags.X GLMFlags.Y GLMFlags.Z GLMFlags.Label]     = textread(TxtName,'%f %f %f %s');
            GLMFlags.Coordinates = 1;
        catch
            [GLMFlags.VOI GLMFlags.Label] = textread(TxtName,'%s %s');
            GLMFlags.Coordinates = 0;
        end
        GLMFlags.project = strrep(strrep(spm_str_manip(TxtName,'rt'),'_VOIs',''),'_ImageList','');
        se_GLM('init')
        
        
    case 'radius'
        GLMFlags.Radius = spm_input('VOI Radius','+0','r',GLMFlags.Radius,1,[2,10]);
        fg = spm_figure('GetWin','Graphics');
        uicontrol(fg,'Style','Pushbutton','Units','normalized','Position',[0.68 .80 .3 .04],'String',{['Sphere: ' num2str(GLMFlags.Radius) 'mm radius']},...
            'HorizontalAlignment','center','FontWeight','bold','Callback','se_GLM(''radius'')','ForegroundColor','k');
        
        
    case 'exit'
        delete(gcf)
        clear
        OverGUI;
end


