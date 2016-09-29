function se_Volumes(varargin)


global NetFlags
clc
if isempty(NetFlags)
    spm_figure('Clear','Graphics');
end

if (nargin==0), Action = 'init'; else, Action = varargin{1}; end

switch lower(Action),
    case 'clears'
        a = zeros(1,numel(NetFlags.Label));
        for i=1:numel(NetFlags.Label)
            a(i) = get(NetFlags.checkBox(i,1),'Value');
        end
        if any(a)
            for i=1:numel(NetFlags.Label)
                set(NetFlags.checkBox(i,1),'Value',0);
            end
        else
            for i=1:numel(NetFlags.Label)
                set(NetFlags.checkBox(i,1),'Value',1);
            end
        end
        
    case 'cleart'
        a = zeros(1,numel(NetFlags.Label));
        for i=1:numel(NetFlags.Label)
            a(i) = get(NetFlags.checkBox(i,2),'Value');
        end
        if any(a)
            for i=1:numel(NetFlags.Label)
                set(NetFlags.checkBox(i,2),'Value',0);
            end
        else
            for i=1:numel(NetFlags.Label)
                set(NetFlags.checkBox(i,2),'Value',1);
            end
        end
        
    case 'init'
        NetFlags.wasison = 2;
        fg = spm_figure('GetWin','Graphics');
        figure(fg);
        if ~isfield(NetFlags,'z')
            try; set(NetFlags.Conf1,'value',1); end
            try; set(NetFlags.Conf2,'value',1); end
            try; set(NetFlags.Conf3,'value',1); end
            try; set(NetFlags.Confx,'value',1); end
        end
        
        Finter = spm('CreateIntWin','on');
        FS   = spm('FontSizes');
        PF   = spm_platform('fonts');
        hFS = FS(get(gcf,'DefaultUicontrolFontSize'));
        a = get(0,'ScreenSize');
        
        if ~isfield(NetFlags,'StatOn');
            NetFlags.StatOn = '0';
        end
        
        ok = 1;
        if ~isfield(NetFlags,'sub')
            NetFlags.Group = 'No group lookup selected';
            ok = 0;
        end
        uicontrol(fg,'Style','Pushbutton','Units','normalized','Position',[0.1 .95 .3 .04],'String',{spm_str_manip(NetFlags.Group,'rt')},...
            'HorizontalAlignment','center','FontWeight','bold','Callback','se_Volumes(''gselect'')','FontSize',10);
        
        if ~isfield(NetFlags,'Label')
            NetFlags.project = 'No VOI file selected';
            NetFlags.Coordinates = 0;
            ok = 0;
        end
        
        uicontrol(fg,'Style','Pushbutton','Units','normalized','Position',[0.1 .9 .3 .04],'String',{spm_str_manip(NetFlags.project)},...
            'HorizontalAlignment','center','FontWeight','bold','Callback','se_Volumes(''vselect'')','FontSize',10);
        
        
        try; sel = get(NetFlags.Mask,'Value'); catch; sel = 1; end
        NetFlags.Mask = uicontrol(fg,'Style','popupmenu','Units','normalized','Position',[0.01 .81 .27 .07],...
            'String',str2mat('VBM 8 Grey Matter','VBM 8 White Matter'),'FontSize',8,'Value',sel,'ToolTipString','Select tissue to analyze');
        
        try; sel = get(NetFlags.mtype,'Value'); catch; sel = 2; end
        NetFlags.mtype = uicontrol(fg,'Style','popupmenu','Units','normalized','Position',[0.28 .81 .39 .07],...
            'String',str2mat('Modulated','Non-linearly modulated','Smoothed Modulated','Smoothed non-linearly modulated'),'FontSize',8,'Value',sel,'ToolTipString','Select type of modulation');
        
        
        try; sel = get(NetFlags.CN,'Value'); catch; sel = 1; end
        NetFlags.CN = uicontrol(fg,'Style','popupmenu','Units','normalized','Position',[0.74 .81 .25 .07],...
            'String',str2mat('No deformation','From ColinMNI','From non-linear MNI152'),'Value',sel,'FontSize',8,'ToolTipString','Move VOIs from standard space');
        
        
        if ok == 0
            uicontrol(fg,'style','pushbutton','units','normalized','position',[0.5 .9 .4 .04],'string',{'Compute/Load volumes'},...
                'horizontalalignment','center','fontweight','bold','callback','se_Volumes(''precomp'')','foregroundcolor','w','FontSize',10);
        else
            uicontrol(fg,'style','pushbutton','units','normalized','position',[0.5 .9 .4 .04],'string',{'Compute/Load volumes'},...
                'horizontalalignment','center','fontweight','bold','callback','se_Volumes(''precomp'')','foregroundcolor','r','FontSize',10);
        end
        
        
        if isfield(NetFlags,'z')
            Label = NetFlags.Label;
            
            
            % Hier jetzt Selection der Confounds fpr raus. max 2? 3?
            try; sel = get(NetFlags.Conf1,'Value'); catch; sel = 1; end
            NetFlags.Conf1 = uicontrol(fg,'Style','popupmenu','Units','normalized','Position',[0.01 .77 .24 .07],...
                'String',{'Confound 1: None' NetFlags.Covariates{:}},'Value',sel,'ToolTipString','Select confound factor',...
                'FontSize',8);
            try; sel = get(NetFlags.Conf2,'Value'); catch; sel = 1; end
            NetFlags.Conf2 = uicontrol(fg,'Style','popupmenu','Units','normalized','Position',[0.25 .77 .24 .07],...
                'String',{'Confound 2: None' NetFlags.Covariates{:}},'Value',sel,'ToolTipString','Select confound factor',...
                'FontSize',8);
            try; sel = get(NetFlags.Conf3,'Value'); catch; sel = 1; end
            NetFlags.Conf3 = uicontrol(fg,'Style','popupmenu','Units','normalized','Position',[0.50 .77 .24 .07],...
                'String',{'Confound 3: None' NetFlags.Covariates{:}},'Value',sel,'ToolTipString','Select confound factor',...
                'FontSize',8);
            try; sel = get(NetFlags.Confx,'Value'); catch; sel = 1; end
            NetFlags.Confx = uicontrol(fg,'Style','popupmenu','Units','normalized','Position',[0.75 .77 .24 .07],...
                'String',{'No removal','1st Order','2nd Order','3rd Order'},'Value',sel,'ToolTipString','Select to remove confounds',...
                'FontSize',8,'foregroundcolor','r','Callback','se_Volumes(''adjust'')');
            
            try; sel = get(NetFlags.Conf1t,'Value'); catch; sel = 1; end
            NetFlags.Conf1t = uicontrol(fg,'Style','popupmenu','Units','normalized','Position',[0.01 .74 .24 .07],...
                'String',str2mat('Categorical','Continuous'),'Value',sel,'ToolTipString','Select confound factor',...
                'FontSize',8);
            try; sel = get(NetFlags.Conf2t,'Value'); catch; sel = 1; end
            NetFlags.Conf2t = uicontrol(fg,'Style','popupmenu','Units','normalized','Position',[0.25 .74 .24 .07],...
                'String',str2mat('Categorical','Continuous'),'Value',sel,'ToolTipString','Select confound factor',...
                'FontSize',8);
            try; sel = get(NetFlags.Conf3t,'Value'); catch; sel = 1; end
            NetFlags.Conf3t = uicontrol(fg,'Style','popupmenu','Units','normalized','Position',[0.50 .74 .24 .07],...
                'String',str2mat('Categorical','Continuous'),'Value',sel,'ToolTipString','Select confound factor',...
                'FontSize',8);

            
            
            
            
            
         uicontrol(fg,'Style','PushButton','Units','normalized','Position',[.02 .02 .08 .05],'ForegroundColor','b','FontWeight','bold',...
            'String','Save','Callback','se_Volumes(''savefiltered'');','ToolTipString','Save current (confound-adjusted) data as Text');
           
          uicontrol(fg,'Style','PushButton','Units','normalized','Position',[.35 .02 .08 .05],'ForegroundColor','g','FontWeight','bold',...
            'String','Clear','Callback','se_Volumes(''clearall'');','ToolTipString','Clear everything');
           
            uicontrol(fg,'Style','Frame','Units','normalized','Position',[.01 .06+.72-.05-.026*numel(Label)-.005 .48 .095+.74-(.06+.72-.024*numel(Label)-.005)   ]);
            uicontrol(fg,'Style','Pushbutton','Units','normalized','Position',[0.02 .05+.74-.05 .46 .04],'ForegroundColor','r',...
                'FontWeight','bold','String','Regions','Callback','se_Volumes(''clearS'')','ToolTipString','Click to (un-) select all');
            
            NetFlags.checkBox = [];
            for i = 1:numel(Label)
                NetFlags.checkBox(i,1) = uicontrol(fg,'Style','checkbox','Units','normalized','Position',[0.05  .06-.05+.723-.026*i .025 .023],'Callback','se_Volumes(''checked'')');
                uicontrol(fg,'Style','text','Units','normalized','Position',[0.1 .06+.72-.05-.026*i .3 .025],'String',{Label{i}},'HorizontalAlignment','left');
                
            end
            
            weiter =  .06+.72-.07-.025*numel(NetFlags.Label);
            
            
            if isfield(NetFlags,'Covariates')
                if numel(NetFlags.Covariates)>0
                    uicontrol(fg,'Style','Pushbutton','Units','normalized','Position',[0.55 .69 .4 .04],'String',{'Covariate analysis'},...
                        'HorizontalAlignment','center','FontWeight','bold','Callback','se_Correlation','ForegroundColor','k');
                end
            end
            
            if isfield(NetFlags,'isPat')
                if numel(unique(NetFlags.isPat))==2
                uicontrol(fg,'Style','Pushbutton','Units','normalized','Position',[0.55 .65 .4 .04],'String',{'Group differences'},...
                    'HorizontalAlignment','center','FontWeight','bold','Callback','se_DualGroup','ForegroundColor','k');
                end
            end
            
            if ~isfield(NetFlags,'pvalue')
                NetFlags.pvalue = 0.05;
            end
            
            uicontrol(fg,'Style','Pushbutton','Units','normalized','Position',[0.55 .55 .4 .04],'String',{['p < ' num2str(NetFlags.pvalue)]},...
                'HorizontalAlignment','center','FontWeight','bold','Callback','se_Volumes(''pvalue'')','ForegroundColor','k');
            
            try; sel = get(NetFlags.Type,'Value'); catch; sel = 1; end
            NetFlags.Type = uicontrol(fg,'Style','popupmenu','Units','normalized','Position',[0.55 .49 .4 .04],'Callback','se_Volumes(''setstat'')',...
                'String',str2mat('T-Tests','Ranksum','Monte-Carlo (Mean)','Monte-Carlo (Median)'),'Value',sel,'ToolTipString','Statistical approach');
            
            try; sel = get(NetFlags.Method,'Value'); catch; sel = 1; end
            NetFlags.Method = uicontrol(fg,'Style','popupmenu','Units','normalized','Position',[0.55 .45 .4 .04],'Callback','se_Volumes(''setstat2'')',...
                'String',str2mat('uncorrected','Bonferroni corrected','FDR corrected'),'Value',sel,'ToolTipString','Correction for multiple comparisons');
   
            try; sel = get(NetFlags.Outlier,'Value'); catch; sel = 1; end
            NetFlags.Outlier = uicontrol(fg,'Style','popupmenu','Units','normalized','Position',[0.55 .41 .4 .04],'Callback','se_Volumes(''adjust'')',...
                'String',str2mat('No outlier exclusion','Exclude FC > 2 std','Exclude FC > 3 std'),'Value',sel,'ToolTipString','Exclude FC outliers for robust regression');
            
            if numel(NetFlags.Label)>1
                uicontrol(fg,'Style','Pushbutton','Units','normalized','Position',[0.55 .37 .4 .04],'String',{'Structural Covariance'},...
                    'HorizontalAlignment','center','FontWeight','bold','Callback','fh_StrucCov(''corr'')','ForegroundColor','k');

                uicontrol(fg,'Style','Pushbutton','Units','normalized','Position',[0.55 .33 .4 .04],'String',{'Hierarchical clustering'},...
                    'HorizontalAlignment','center','FontWeight','bold','Callback','fh_StrucCov(''ward'')','ForegroundColor','k');
                
                try; sel = get(NetFlags.Corr,'Value'); catch; sel = 1; end
                NetFlags.Corr = uicontrol(fg,'Style','popupmenu','Units','normalized','Position',[0.55 .27 .4 .04],'Callback','se_Volumes(''setstat'')',...
                    'String',str2mat('Pearson linear correlation','Spearman rank correlation','Partial correlation'),'Value',sel,'ToolTipString','Type of Correlation');

            end
            
        end
        
        uicontrol(fg,'Style','PushButton','Units','normalized','Position',[.9 .02 .08 .05],'ForegroundColor','r','FontWeight','bold',...
            'String','EXIT','Callback','se_Volumes(''exit'');','ToolTipString','quit');
        F = findall(allchild(0),'Flat','Tag','Show Effect');
        close(F)
        F = findall(allchild(0),'Flat','Tag','Results');
        
        close(F)
        
        
        
    case 'adjust'
        
        NetFlags.z = NetFlags.rawz;
        NetFlags.Cov = NetFlags.rawCov;        

        if get(NetFlags.Confx,'Value')>1
            
        %%% mean center outlieres for confound regression
            NoutC = repmat({''},numel(NetFlags.Label));
            NoutP = repmat({''},numel(NetFlags.Label));
            if get(NetFlags.Outlier,'Value')>=2
                for i=1:numel(NetFlags.Label)
                    NoutC{i} = find(abs(NetFlags.z(NetFlags.isPat==0,i))>mean(NetFlags.z(NetFlags.isPat==0,i))+get(NetFlags.Outlier,'Value')*std(NetFlags.z(NetFlags.isPat==0,i)));
                    NoutP{i} = find(abs(NetFlags.z(NetFlags.isPat==1,i))>mean(NetFlags.z(NetFlags.isPat==1,i))+get(NetFlags.Outlier,'Value')*std(NetFlags.z(NetFlags.isPat==1,i)));
                    NetFlags.z(NoutC{i},i) = mean(NetFlags.z(NetFlags.isPat==0,i));
                    NetFlags.z(NoutP{i},i) = mean(NetFlags.z(NetFlags.isPat==1,i));
                end
                cnt = 0;
                for i=1:numel(NetFlags.Cov(1,:))
                    CoutC{i} = find(abs(NetFlags.Cov(NetFlags.isPat==0,i))>mean(NetFlags.Cov(NetFlags.isPat==0,i))+get(NetFlags.Outlier,'Value')*std(NetFlags.Cov(NetFlags.isPat==0,i)));
                    CoutP{i} = find(abs(NetFlags.Cov(NetFlags.isPat==1,i))>mean(NetFlags.Cov(NetFlags.isPat==1,i))+get(NetFlags.Outlier,'Value')*std(NetFlags.Cov(NetFlags.isPat==1,i)));
                    NetFlags.Cov(CoutC{i},i) = mean(NetFlags.Cov(NetFlags.isPat==0,i));
                    NetFlags.Cov(CoutP{i},i) = mean(NetFlags.Cov(NetFlags.isPat==1,i));
                    cnt = cnt + numel(CoutC{i}>0);
                    cnt = cnt + numel(CoutP{i}>0);
                end
            else
                CoutC = repmat({''},numel(NetFlags.Cov(1,:)));
                CoutP = repmat({''},numel(NetFlags.Cov(1,:)));
                for i=1:numel(NetFlags.Cov(1,:))
                    CoutC{i} = '';
                    CoutP{i} = '';
                end
            end
            
            Conf = [get(NetFlags.Conf1,'Value') get(NetFlags.Conf2,'Value') get(NetFlags.Conf3,'Value')]; 
            Conft = [get(NetFlags.Conf1t,'Value') get(NetFlags.Conf2t,'Value') get(NetFlags.Conf3t,'Value')]; 
            
            [Conf ia] = unique(Conf(Conf>1));
            Conf = Conf-1;
            cont = find((Conft(ia)-1)>0);

            modl = get(NetFlags.Confx,'Value')-1;
            
            for i=1:size(NetFlags.z,2)
                fprintf(1,'%s\n',[int2str(i) '/' int2str(size(NetFlags.z,2)) ': ' NetFlags.Label{i}])
                c = mean(NetFlags.z(:,i));
                [p,table,stats] = anovan(NetFlags.z(:,i)-c,mat2cell(NetFlags.Cov(:,Conf),size(NetFlags.Cov,1),[ones(1,numel(Conf))]),...
                    'model',modl,'display','off','sstype',3,'continuous',cont);
                NetFlags.z(:,i) = stats.resid+c;
            end
            
            for i=1:size(NetFlags.Cov,2)
                if ~any(Conf==i)
                    fprintf(1,'%s\n',[int2str(i) '/' int2str(size(NetFlags.Cov,2)) ': ' NetFlags.Covariates{i}])
                    Q = ~isnan(NetFlags.Cov(:,i));
                    c = mean(NetFlags.Cov(Q,i));
                    [p,table,stats] = anovan(NetFlags.Cov(Q,i)-c,mat2cell(NetFlags.Cov(Q,Conf),sum(Q),[ones(1,numel(Conf))]),...
                        'model',get(NetFlags.Confx,'Value')-1,'display','off','sstype',3,'continuous',cont);
                    NetFlags.Cov(Q,i) = stats.resid+c;
                end
            end
            
            for i=1:numel(NetFlags.Label)
                NetFlags.z(NoutC{i},i)=NaN;
                NetFlags.z(NoutP{i},i)=NaN;
                NetFlags.Cov(CoutC{i},i) = NaN;
                NetFlags.Cov(CoutP{i},i) = NaN;
            end
            if get(NetFlags.Outlier,'Value')>=2
                spm('alert"',{['Confounds removed with ' int2str((sum(isnan(NetFlags.z(:)))/numel(NetFlags.z))*100) ' % Volumes & ' ...
                    int2str((cnt/numel(NetFlags.Cov))*100) ' % Covariates outlieres excluded']});
            else
            spm('alert"',{'Confounds removed without outlier exclusion!'});
            end            
        else
            if get(NetFlags.Outlier,'Value')>=2
                for i=1:numel(NetFlags.Label)
                   NetFlags.z(find(abs(NetFlags.z(NetFlags.isPat==0,i))>mean(NetFlags.z(NetFlags.isPat==0,i))+get(NetFlags.Outlier,'Value')*std(NetFlags.z(NetFlags.isPat==0,i))),i) = NaN;
                   NetFlags.z(find(abs(NetFlags.z(NetFlags.isPat==1,i))>mean(NetFlags.z(NetFlags.isPat==1,i))+get(NetFlags.Outlier,'Value')*std(NetFlags.z(NetFlags.isPat==1,i))),i) = NaN;
                end
                cnt = 0;
                for i=1:numel(NetFlags.Cov(1,:))
                    CoutC{i} = find(abs(NetFlags.Cov(NetFlags.isPat==0,i))>mean(NetFlags.Cov(NetFlags.isPat==0,i))+get(NetFlags.Outlier,'Value')*std(NetFlags.Cov(NetFlags.isPat==0,i)));
                    CoutP{i} = find(abs(NetFlags.Cov(NetFlags.isPat==1,i))>mean(NetFlags.Cov(NetFlags.isPat==1,i))+get(NetFlags.Outlier,'Value')*std(NetFlags.Cov(NetFlags.isPat==1,i)));
                    NetFlags.Cov(CoutC{i},i) = NaN;
                    NetFlags.Cov(CoutP{i},i) = NaN;
                    cnt = cnt + numel(CoutC{i}>0);
                    cnt = cnt + numel(CoutP{i}>0);
                end
                spm('alert"',{[int2str((sum(isnan(NetFlags.z(:)))/numel(NetFlags.z))*100) ' % Volumes & ' ...
                    int2str((cnt/numel(NetFlags.Cov))*100) ' % Covariates outlieres excluded']});
            else
%             spm('alert"',{'No outliers removed!'});                
            end            
        end
        
        switch NetFlags.StatOn
            case 'S'
                se_SingleGroup
            case 'G'
                se_DualGroup
            case 'C'
                se_Correlation
        end
        
        
        
        
    case 'significance'
        if get(NetFlags.Type,'value')==3
            spm('alert*',{'No Monto-Carlo for baseline tests'})
        else
            se_SingleGroup
        end
        
        
        
    case 'setstat'
        F = findall(allchild(0),'Flat','Tag','Show Effect');
        close(F)
        switch NetFlags.StatOn
            case 'S'
                if get(NetFlags.Type,'value')<3;  se_SingleGroup; end
            case 'G'
                se_DualGroup
            case 'C'
                se_Correlation
        end
        
        
    case 'setstat2'
        F = findall(allchild(0),'Flat','Tag','Show Effect');
        close(F)
        switch NetFlags.StatOn
            case 'S'
                se_SingleGroup
            case 'G'
                se_DualGroup
            case 'C'
                se_Correlation
        end
        
        
    case 'pvalue'
        NetFlags.pvalue = spm_input('p value','+0','r',NetFlags.pvalue,1,[0,1]);
        fg = spm_figure('GetWin','Graphics');
        uicontrol(fg,'Style','Pushbutton','Units','normalized','Position',[0.55 .55 .4 .04],'String',{['p < ' num2str(NetFlags.pvalue)]},...
            'HorizontalAlignment','center','FontWeight','bold','Callback','se_Volumes(''pvalue'')','ForegroundColor','k');
        F = findall(allchild(0),'Flat','Tag','Show Effect');
        close(F)
        switch NetFlags.StatOn
            case 'S'
                se_SingleGroup
            case 'G'
                se_DualGroup
            case 'C'
                se_Correlation
        end
        
        
    case 'precomp'
        ok = [0 0];
        if isfield(NetFlags,'sub')
            ok(1) = 1;
        else
            spm('alert*',{'No group lookup selected'})
        end
        if isfield(NetFlags,'Label')
            ok(2) = 1;
        else
            spm('alert*',{'No VOI file selected'})
        end
        if sum(ok)==2
            GM   = {'VBM8_GM','VBM8_WM'};
            CN = {'Native','Colin2VBM'};
            MT = {'m','m0','sm','sm0'};
            des = [MT{get(NetFlags.mtype,'value')} '_' GM{get(NetFlags.Mask,'value')} '_' CN{get(NetFlags.CN,'value')}];
            
            try
                tmp = load(fullfile(pwd,'Projects','VBM','Regional',spm_str_manip(NetFlags.Group,'rt'),[NetFlags.project '_' des '.mat']));
                NetFlags.z = tmp.z; clear tmp
                
            catch
                se_precomputeVol(des)
            end
            
            if numel(NetFlags.Label)==2
                if strcmpi(questdlg('Include Differntial ?','Interaction','Yes','No','No'),'yes')
                    NetFlags.Label{3} = [NetFlags.Label{1} ' - ' NetFlags.Label{2}];
                    NetFlags.z(:,3) = NetFlags.z(:,1) - NetFlags.z(:,2);
                end
            end
            
            NetFlags.rawz = NetFlags.z;
            se_Volumes('init')
            se_Volumes('adjust');
            
        end
        
        
        
    case 'mprecomp'
        mgroup = spm_select(Inf,'any',['Select group descriptions'],[],fullfile(pwd,'Lookups'),'.mat',1);
        mvoi   = spm_select(Inf,'any',['Select coordinate list'],[],fullfile(pwd,'VOIs'),'.txt',1);
        
        for gg = 1:size(mgroup,1)
            for vv = 1:size(mvoi,1)
                
                
                NetFlags.Group = deblank(mgroup(gg,:));
                tmp = load(NetFlags.Group);
                NetFlags.Covariates = tmp.Covariates;
                NetFlags.sub = tmp.sub(tmp.hasVBM>0);
                NetFlags.Cov = tmp.Cov(tmp.hasVBM>0,:);
                NetFlags.isPat = tmp.isPat(tmp.hasVBM>0);
                NetFlags.SubDir = tmp.SubDir(tmp.hasVBM>0);
                
                clear tmp
                
                TxtName = deblank(mvoi(vv,:));
                try
                    [NetFlags.VOI NetFlags.Label] = textread(TxtName,'%s %s');
                    NetFlags.Coordinates = 0;
                catch
                    spm('alert*',{'Not a valid image-VOI file'})
                    error
                end
                NetFlags.project = strrep(strrep(spm_str_manip(TxtName,'rt'),'_VOIs',''),'_ImageList','');
                NetFlags.Nodes = zeros(numel(NetFlags.Label,2));
                
                
                for gmm = 1:2
                        for ms = 1:4
                            set(NetFlags.Mask,'value',gmm); set(NetFlags.mtype,'value',ms);
                            GM   = {'GM','WM'};
                            Glob = {'TTV','TBV','NN'};
                            CN = {'Native','Colin2VBM'};
                            MT = {'m','m0','sm','sm0'};
                            des = [MT{get(NetFlags.mtype,'value')} '_' GM{get(NetFlags.Mask,'value')} '_' CN{get(NetFlags.CN,'value')} ];
                            
                            se_precomputeVol(des)
                        end
                end
                
                
                
                
            end
        end
        
        
    case 'gselect'
        NetFlags.Group = spm_select(1,'mat',['Select group description'],[],fullfile(pwd,'Lookups'),'.mat',1);
        tmp = load(NetFlags.Group);
        NetFlags.Covariates = tmp.Covariates;
        NetFlags.sub = tmp.sub(tmp.hasVBM==1);
        try; NetFlags.Cov = tmp.Cov(tmp.hasVBM==1,:); end
        try; NetFlags.rawCov = tmp.Cov(tmp.hasVBM==1,:); end
        NetFlags.isPat = tmp.isPat(tmp.hasVBM==1);
        NetFlags.isPat(isnan(NetFlags.isPat)) = 0;
        NetFlags.SubDir = tmp.SubDir(tmp.hasVBM==1);
        clear tmp
        if isfield(NetFlags,'z'); NetFlags = rmfield(NetFlags,'z'); end
        try; NetFlags = rmfield(NetFlags,'activeCov'); end
        try; NetFlags = rmfield(NetFlags,'StatOn');end
        try; NetFlags = rmfield(NetFlags,'Connections');end
        try; NetFlags = rmfield(NetFlags,'Results');end
        try; NetFlags = rmfield(NetFlags,'Description');end
        try; NetFlags = rmfield(NetFlags,'Type');end
        
        try, NetFlags = rmfield(NetFlags,'z'); end
        se_Volumes('init')
        se_Volumes('vselect')
        
        
    case 'vselect'
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
        
        
        NetFlags.Nodes = zeros(numel(NetFlags.Label,2));
        if isfield(NetFlags,'z'); NetFlags = rmfield(NetFlags,'z'); end
        try; NetFlags = rmfield(NetFlags,'activeCov'); end
        try; NetFlags = rmfield(NetFlags,'StatOn');end
        try; NetFlags = rmfield(NetFlags,'Connections');end
        try; NetFlags = rmfield(NetFlags,'Results');end
        try; NetFlags = rmfield(NetFlags,'Description');end
        try; NetFlags = rmfield(NetFlags,'Type');end
        try, NetFlags = rmfield(NetFlags,'z'); end
        se_Volumes('init')
        
        
    case 'exit'
        delete(gcf)
%         se_figure('Clear','Graphics');
        clear
        OverGUI;
        
        
        
        
    case 'savefiltered'
        
        [file, pfad] = uiputfile({'*.txt','Textfile (*.txt)'},'Save adjusted data as text in',fullfile(pwd,'Projects','VBM','Regional',spm_str_manip(NetFlags.Group,'rt')));
        
        wrz = [];
        wn  = {};
        cnt = 1;
        
        for i=1:numel(NetFlags.Label)
            wrz(cnt,:) = NetFlags.z(:,i)';
            wn{cnt} = [NetFlags.Label{i}];
            cnt = cnt+1;
        end
        wrz = wrz';
        
        fid = fopen(fullfile(pfad,file),'w+');
        fprintf(fid,'%s\t','Subject');
        for i=1:numel(wn)
            fprintf(fid,'%s\t',wn{i});
        end
        for ii=1:numel(NetFlags.sub)
            fprintf(fid,'\n%s\t',NetFlags.sub{ii});
            
            for i=1:numel(wn)
                fprintf(fid,'%8.6f\t',wrz(ii,i));
            end
        end
        fclose(fid);
        
        
        fid = fopen(fullfile(pfad,strrep(file,'.','_sub.')),'w+');
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
        

        
    case 'clearall'
        clear global
%        se_figure('Clear','Graphics');
        se_Volumes('init')
       
end


