function se_Network(varargin)


global NetFlags


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
        NetFlags.wasison = 1;
        try
            NetFlags.Mask = get(NetFlags.Mask,'Value');
            NetFlags.PCA  = get(NetFlags.PCA,'Value');
            NetFlags.Glob = get(NetFlags.Glob,'Value');
            NetFlags.Filt = get(NetFlags.Filt,'Value');
        end
        fg = spm_figure('GetWin','Graphics');
        if ~isfield(NetFlags,'z')
            figure(fg); clf;
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
            'HorizontalAlignment','center','FontWeight','bold','Callback','se_Network(''gselect'')','FontSize',10,...
            'ToolTipString','Choose Lookup-file describing resting-state sample');
        
        if ~isfield(NetFlags,'Label')
            NetFlags.project = 'No VOI file selected';
            NetFlags.Coordinates = 0;
            ok = 0;
        end
        
        uicontrol(fg,'Style','Pushbutton','Units','normalized','Position',[0.1 .9 .3 .04],'String',{spm_str_manip(NetFlags.project)},...
            'HorizontalAlignment','center','FontWeight','bold','Callback','se_Network(''vselect'')','FontSize',10,...
            'ToolTipString','File may contain coordinates or image VOIs');
        
        
        try; sel = NetFlags.Mask; catch; sel = 1; end
        NetFlags.Mask = uicontrol(fg,'Style','popupmenu','Units','normalized','Position',[0.02 .8 .18 .08],...
            'String',str2mat('VOIs: No GM mask','Individual median-split','Group median-split'),'FontSize',8,'Value',sel,...
            'ToolTipString','Select method for VOI gray matter masking');
        
        try; if mean(NetFlags.fix)>.7
                NetFlags.PCA = uicontrol(fg,'Style','popupmenu','Units','normalized','Position',[0.21 .8 .18 .08],...
                    'String',str2mat('PCA denoising','No PCA denoising','FIX de-noising'),'Value',NetFlags.PCA,'foregroundcolor','r','FontSize',10,...
                    'ToolTipString','FMRIB''s ICA-based X-noiseifier');
            else
                NetFlags.PCA = uicontrol(fg,'Style','popupmenu','Units','normalized','Position',[0.21 .8 .18 .08],...
                    'String',str2mat('PCA denoising','No PCA denoising'),'Value',NetFlags.PCA,'FontSize',8,'ToolTipString','Include 5 component PCA denoising');
            end
        catch
            try; sel = NetFlags.PCA; catch; sel = 2; end
            NetFlags.PCA = uicontrol(fg,'Style','popupmenu','Units','normalized','Position',[0.21 .8 .18 .08],...
                'String',str2mat('PCA denoising','No PCA denoising'),'Value',sel,'FontSize',8,'ToolTipString','Include 5 component PCA denoising');
        end
        
        try; sel = NetFlags.Glob; catch; try; if any(NetFlags.isPat)>0; sel = 7; else sel = 6; end; catch; sel = 7; end; end
        NetFlags.Glob = uicontrol(fg,'Style','popupmenu','Units','normalized','Position',[0.4 .8 .22 .08],'String',str2mat...
            ('Full tissue signal regression','Linear tissue signal regression','Full WM/CSF signal regression','Linear WM/CSF signal regression',...
            'Full global signal regression','Linear global signal regression', 'No Global signal regression'),...
            'Value',sel,'FontSize',8,'ToolTipString','Mean signal regression');
        
        try; sel = NetFlags.Filt; catch; sel = 1; end
        NetFlags.Filt = uicontrol(fg,'Style','popupmenu','Units','normalized','Position',[0.63 .8 .13 .08],...
            'String',str2mat('Band-pass filter','High-pass filter'),'Value',sel,'FontSize',8,'ToolTipString','Filter for BOLD frequencies');
        
        
        if NetFlags.Coordinates>0
            if ~isfield(NetFlags,'Radius'); NetFlags.Radius = 5; end
            uicontrol(fg,'Style','Pushbutton','Units','normalized','Position',[0.78 .84 .2 .04],'String',{['Sphere: ' num2str(NetFlags.Radius) 'mm radius']},...
                'HorizontalAlignment','center','FontWeight','bold','Callback','se_Network(''radius'')','ForegroundColor','k','FontSize',8);
        end
        
        
        
        
        
        if ok == 0
            uicontrol(fg,'style','pushbutton','units','normalized','position',[0.5 .9 .4 .04],'string',{'Compute/Load connectivities'},...
                'callback','se_Network(''precomp'')','foregroundcolor','w','FontSize',10);
        else
            uicontrol(fg,'style','pushbutton','units','normalized','position',[0.5 .9 .4 .04],'string',{'Compute/Load connectivities'},...
                'horizontalalignment','center','fontweight','bold','callback','se_Network(''precomp'')','foregroundcolor','r','FontSize',10);
        end
        
        uicontrol(fg,'Style','Pushbutton','Units','normalized','Position',[0.5 .95 .4 .04],'String',{'Pre-compute all variants'},...
            'HorizontalAlignment','center','FontWeight','bold','Callback','se_Network(''mprecomp'')','FontSize',10);
        
        if isfield(NetFlags,'z')
            Label = NetFlags.Label;
           
            
          
            try; sel = get(NetFlags.CORR,'Value'); catch; sel = 1; end
            NetFlags.CORR = uicontrol(fg,'Style','popupmenu','Units','normalized','Position',[0.25 .805 .44 .05],'String',...
                str2mat('                            Pearson correlation','                     Spearmann rank correlation'...
                    ,'                             Partial correlation','                    Regularised partial correlation'),...
                'fontweight','bold','Value',sel,'FontSize',8,'ToolTipString','Type of correlation','foregroundcolor','r','Callback','se_Network(''adjust'')');

                       
            
            % Hier jetzt Selection der Confounds fpr raus. max 2? 3?
            try; sel = get(NetFlags.Conf1,'Value'); catch; sel = 1; end
            NetFlags.Conf1 = uicontrol(fg,'Style','popupmenu','Units','normalized','Position',[0.01 .76 .24 .07],...
                'String',{'Confound 1: None' NetFlags.Covariates{:}},'Value',sel,'ToolTipString','Select confound factor',...
                'FontSize',8);
            try; sel = get(NetFlags.Conf2,'Value'); catch; sel = 1; end
            NetFlags.Conf2 = uicontrol(fg,'Style','popupmenu','Units','normalized','Position',[0.26 .76 .24 .07],...
                'String',{'Confound 2: None' NetFlags.Covariates{:}},'Value',sel,'ToolTipString','Select confound factor',...
                'FontSize',8);
            try; sel = get(NetFlags.Conf3,'Value'); catch; sel = 1; end
            NetFlags.Conf3 = uicontrol(fg,'Style','popupmenu','Units','normalized','Position',[0.51 .76 .24 .07],...
                'String',{'Confound 3: None' NetFlags.Covariates{:}},'Value',sel,'ToolTipString','Select confound factor',...
                'FontSize',8);

            try; sel = get(NetFlags.Conf1t,'Value'); catch; sel = 1; end
            NetFlags.Conf1t = uicontrol(fg,'Style','popupmenu','Units','normalized','Position',[0.01 .74 .24 .07],...
                'String',str2mat('Categorical','Continuous'),'Value',sel,'ToolTipString','Select confound factor',...
                'FontSize',8);
            try; sel = get(NetFlags.Conf2t,'Value'); catch; sel = 1; end
            NetFlags.Conf2t = uicontrol(fg,'Style','popupmenu','Units','normalized','Position',[0.26 .74 .24 .07],...
                'String',str2mat('Categorical','Continuous'),'Value',sel,'ToolTipString','Select confound factor',...
                'FontSize',8);
            try; sel = get(NetFlags.Conf3t,'Value'); catch; sel = 1; end
            NetFlags.Conf3t = uicontrol(fg,'Style','popupmenu','Units','normalized','Position',[0.51 .74 .24 .07],...
                'String',str2mat('Categorical','Continuous'),'Value',sel,'ToolTipString','Select confound factor',...
                'FontSize',8);


            try; sel = get(NetFlags.Confx,'Value'); catch; sel = 1; end
            NetFlags.Confx = uicontrol(fg,'Style','popupmenu','Units','normalized','Position',[0.77 .76 .2 .07],...
                'String',{'No removal','1st Order','2nd Order','3rd Order'},'Value',sel,'ToolTipString','Select to remove confounds',...
                'FontSize',9,'foregroundcolor','r','Callback','se_Network(''adjust'')');
            
            
             
         uicontrol(fg,'Style','PushButton','Units','normalized','Position',[.02 .00 .08 .05],'ForegroundColor','b','FontWeight','bold',...
            'String','Save','Callback','se_Network(''savefiltered'');','ToolTipString','Save current (confound-adjusted) data as Text');
           
         uicontrol(fg,'Style','PushButton','Units','normalized','Position',[.35 .00 .08 .05],'ForegroundColor','g','FontWeight','bold',...
            'String','Clear','Callback','se_Network(''clearall'');','ToolTipString','Clear everything');

        
        
            
            uicontrol(fg,'Style','Frame','Units','normalized','Position',[.01 .06+.72-.05-.024*numel(Label)-.005 .48 .095+.74-(.06+.72-.024*numel(Label)-.005)   ]);
            uicontrol(fg,'Style','Pushbutton','Units','normalized','Position',[0.02 .05+.74-.05 .46 .04],'ForegroundColor','r',...
                'FontWeight','bold','String','Seeds','Callback','se_Network(''clearS'')','ToolTipString','Click to (un-) select all');
            
            uicontrol(fg,'Style','Frame','Units','normalized','Position',[.51 .06-.05+.72-.024*numel(Label)-.005 .48 .095+.74-(.06+.72-.024*numel(Label)-.005)   ]);
            uicontrol(fg,'Style','Pushbutton','Units','normalized','Position',[0.52  .05+.74-.05 .46 .04],'ForegroundColor','r',...
                'FontWeight','bold','String','Targets','Callback','se_Network(''clearT'')','ToolTipString','Click to (un-) select all');
            
            for i = 1:numel(Label)
                NetFlags.checkBox(i,1) = uicontrol(fg,'Style','checkbox','Units','normalized','Position',[0.05  .06-.05+.723-.024*i .025 .023],'Callback','se_Network(''checked'')');
                uicontrol(fg,'Style','text','Units','normalized','Position',[0.1 .06+.72-.05-.024*i .3 .023],'String',{Label{i}},'HorizontalAlignment','left');
                
                NetFlags.checkBox(i,2) =   uicontrol(fg,'Style','checkbox','Units','normalized','Position',[0.55 .06-.05+.723-.024*i .025 .023],'Callback','se_Network(''checked'')');
                uicontrol(fg,'Style','text','Units','normalized','Position',[0.6  .06+.72-.05-.024*i .3 .023],'String',{Label{i}},'HorizontalAlignment','left');
                
            end
            
            weiter =  .06+.72-.07-.025*numel(NetFlags.Label);
            
            
            
            if any(NetFlags.isPat)>0
                try; sel = get(NetFlags.Sing,'Value'); catch; sel = 1; end
                NetFlags.Sing = uicontrol(fg,'Style','popupmenu','Units','normalized','Position',[0.05 weiter-.03 .4 .04],'Callback','se_Network(''significance'')',...
                    'String',str2mat('Test significant connections','Only Controls','Only Patients'),'Value',sel,'FontWeight','bold');
            else
                uicontrol(fg,'Style','Pushbutton','Units','normalized','Position',[0.05 weiter-.03 .4 .04],'String',{'Test significant connections'},...
                    'HorizontalAlignment','center','FontWeight','bold','Callback','se_Network(''significance'')','ForegroundColor','k');
                NetFlags.Sing = 1;
            end
            
            uicontrol(fg,'Style','Pushbutton','Units','normalized','Position',[0.05 weiter-.08 .4 .04],'String',{'Hierarchical clustering'},...
                'HorizontalAlignment','center','FontWeight','bold','ForegroundColor','k','Callback','fh_FSLNets','ToolTipString','WARD clustering over all VOIs');

            
            if any(NetFlags.isPat)>0
                uicontrol(fg,'Style','Pushbutton','Units','normalized','Position',[0.05 weiter-.13 .4 .04],'String',{'Group differences'},...
                    'HorizontalAlignment','center','FontWeight','bold','Callback','se_DualGroup','ForegroundColor','k');
                weiter = weiter-.05;
            end
            
                        
            if isfield(NetFlags,'Covariates')
                if numel(NetFlags.Covariates)>0
                    uicontrol(fg,'Style','Pushbutton','Units','normalized','Position',[0.05 weiter-.13 .4 .04],'String',{'Covariate analysis'},...
                        'HorizontalAlignment','center','FontWeight','bold','Callback','se_Correlation','ForegroundColor','k');
                end
            end
            


            
            if ~isfield(NetFlags,'pvalue')
                NetFlags.pvalue = 0.05;
            end
            
            uicontrol(fg,'Style','Pushbutton','Units','normalized','Position',[0.55 weiter-.03 .4 .04],'String',{['p < ' num2str(NetFlags.pvalue)]},...
                'HorizontalAlignment','center','FontWeight','bold','Callback','se_Network(''pvalue'')','ForegroundColor','k');
            
            try; sel = get(NetFlags.Type,'Value'); catch; sel = 1; end
            NetFlags.Type = uicontrol(fg,'Style','popupmenu','Units','normalized','Position',[0.55 weiter-.08 .4 .04],'Callback','se_Network(''setstat'')',...
                'String',str2mat('T-Tests','Ranksum','Monte-Carlo (Mean)','Monte-Carlo (Median)'),'Value',sel,'ToolTipString','Statistical approach');
            
            try; sel = get(NetFlags.Method,'Value'); catch; sel = 1; end
            NetFlags.Method = uicontrol(fg,'Style','popupmenu','Units','normalized','Position',[0.55 weiter-.11 .4 .04],'Callback','se_Network(''setstat2'')',...
                'String',str2mat('uncorrected','Bonferroni corrected','FDR corrected'),'Value',sel,'ToolTipString','Correction for multiple comparisons');
            
            try; sel = get(NetFlags.Outlier,'Value'); catch; sel = 1; end
            NetFlags.Outlier = uicontrol(fg,'Style','popupmenu','Units','normalized','Position',[0.55 weiter-.147 .4 .04],'Callback','se_Network(''adjust'')',...
                'String',str2mat('No outlier exclusion','Exclude FC > 2 std','Exclude FC > 3 std'),'Value',sel,'ToolTipString','Exclude FC outliers for robust regression');
            
        end
        
        uicontrol(fg,'Style','PushButton','Units','normalized','Position',[.9 .00 .08 .05],'ForegroundColor','r','FontWeight','bold',...
            'String','EXIT','Callback','se_Network(''exit'');','ToolTipString','quit');
        F = findall(allchild(0),'Flat','Tag','Show Effect');
        close(F)
        F = findall(allchild(0),'Flat','Tag','Results');
        
        close(F)
        
    case 'adjust'
        
        if get(NetFlags.CORR,'Value')==1
                NetFlags.z = NetFlags.origCorr;
        elseif get(NetFlags.CORR,'Value')==2
                NetFlags.z = NetFlags.origRanC;
        elseif get(NetFlags.CORR,'Value')==3
                NetFlags.z = NetFlags.origParC;
        elseif get(NetFlags.CORR,'Value')==4
                NetFlags.z = NetFlags.origRegPC;
        end
        
        NetFlags.Cov = NetFlags.rawCov;
        
        tmp = load(NetFlags.Group);
        
        if get(NetFlags.PCA,'value')==3
            incl=find(tmp.hasRSfix==1);
        else
            incl=find(tmp.hasRS==1);
        end


        NetFlags.sub = tmp.sub(incl);
        NetFlags.Cov = tmp.Cov(incl,:);
        NetFlags.rawCov = tmp.Cov(incl,:);
        NetFlags.isPat = tmp.isPat(incl);
        NetFlags.SubDir = tmp.SubDir(incl);
        try NetFlags.fix = tmp.hasRSfix(incl); end
        try NetFlags.matchQ=tmp.matchQ(incl); end
        
        if get(NetFlags.Confx,'Value')>1
            
        %%% mean center outlieres for confound regression
            NoutC = repmat({''},numel(NetFlags.Label));
            NoutP = repmat({''},numel(NetFlags.Label));
            if get(NetFlags.Outlier,'Value')>=2
                for i=1:numel(NetFlags.Label)
                    for ii=1:numel(NetFlags.Label)
                       meanC = mean(NetFlags.z(NetFlags.isPat==0,i,ii));
                       stdC  = std(NetFlags.z(NetFlags.isPat==0,i,ii));
                       meanP = mean(NetFlags.z(NetFlags.isPat==1,i,ii));
                       stdP  = std(NetFlags.z(NetFlags.isPat==1,i,ii));
                       NoutC{i,ii} = find((NetFlags.z(:,i,ii) > meanC+get(NetFlags.Outlier,'Value')*stdC & NetFlags.isPat(:)==0) | ...
                                          (NetFlags.z(:,i,ii) < meanC-get(NetFlags.Outlier,'Value')*stdC & NetFlags.isPat(:)==0));
                       NoutP{i,ii} = find((NetFlags.z(:,i,ii) > meanP+get(NetFlags.Outlier,'Value')*stdP & NetFlags.isPat(:)==1) | ...
                                          (NetFlags.z(:,i,ii) < meanP-get(NetFlags.Outlier,'Value')*stdP & NetFlags.isPat(:)==1));
                       NetFlags.z(NoutC{i,ii},i,ii) = mean(NetFlags.z(NetFlags.isPat==0,i,ii));
                       NetFlags.z(NoutP{i,ii},i,ii) = mean(NetFlags.z(NetFlags.isPat==1,i,ii));
                    end
                end
            end

        %%% Do confound regression
            Conf = [get(NetFlags.Conf1,'Value') get(NetFlags.Conf2,'Value') get(NetFlags.Conf3,'Value')]; 
            Conft = [get(NetFlags.Conf1t,'Value') get(NetFlags.Conf2t,'Value') get(NetFlags.Conf3t,'Value')]; 
            
            [Conf ia] = unique(Conf(Conf>1));
            Conf = Conf-1;
            cont = (Conft(ia)-1)>0;
            
            for i=1:size(NetFlags.z,2)-1
                for ii=i+1:size(NetFlags.z,2)
                    c = mean(NetFlags.z(:,i,ii));
                    [p,table,stats] = anovan(NetFlags.z(:,i,ii)-c,mat2cell(NetFlags.Cov(:,Conf),size(NetFlags.Cov,1),[ones(1,numel(Conf))]),...
                        'model',get(NetFlags.Confx,'Value')-1,'display','off','sstype',3,'continuous',cont);
                    NetFlags.z(:,i,ii) = stats.resid+c;
                    NetFlags.z(:,ii,i) = stats.resid+c;
                    
                end
            end
            
            for i=1:size(NetFlags.rawCov,2)
                if ~any(Conf==i)
                    Q = ~isnan(NetFlags.rawCov(:,i));
                    c = mean(NetFlags.rawCov(Q,i));
                    [p,table,stats] = anovan(NetFlags.rawCov(Q,i)-c,mat2cell(NetFlags.Cov(Q,Conf),sum(Q),[ones(1,numel(Conf))]),...
                        'model',get(NetFlags.Confx,'Value')-1,'display','off','sstype',3,'continuous',cont);
                    NetFlags.Cov(Q,i) = stats.resid+c;
                else
                    NetFlags.Cov(:,i) = NetFlags.rawCov(:,i);
                end
            end
            
        %%% Exclude outlieres
            for i=1:numel(NetFlags.Label)
                for ii=1:numel(NetFlags.Label)
                    NetFlags.z(NoutC{i,ii},i,ii)=NaN;
                    NetFlags.z(NoutP{i,ii},i,ii)=NaN;
                end
            end
            
            if get(NetFlags.Outlier,'Value')>=2
                spm('alert"',{['Confounds removed with ' int2str((sum(isnan(NetFlags.z(:)))/numel(NetFlags.z))*100) ' % FC outlieres excluded']});
            else
            spm('alert"',{'Confounds removed without outlier exclusion!'});
            end
        else
            
        %%% Exclude outlieres
            if get(NetFlags.Outlier,'Value')>=2
                for i=1:numel(NetFlags.Label)
                    for ii=1:numel(NetFlags.Label)
                       meanC = mean(NetFlags.z(NetFlags.isPat==0,i,ii));
                       stdC  = std(NetFlags.z(NetFlags.isPat==0,i,ii));
                       meanP = mean(NetFlags.z(NetFlags.isPat==1,i,ii));
                       stdP  = std(NetFlags.z(NetFlags.isPat==1,i,ii));
                       NetFlags.z(find(NetFlags.z(:,i,ii) > meanC+get(NetFlags.Outlier,'Value')*stdC & NetFlags.isPat(:)==0) ,i,ii) = NaN;
                       NetFlags.z(find(NetFlags.z(:,i,ii) < meanC-get(NetFlags.Outlier,'Value')*stdC & NetFlags.isPat(:)==0) ,i,ii) = NaN;
                       NetFlags.z(find(NetFlags.z(:,i,ii) > meanP+get(NetFlags.Outlier,'Value')*stdP & NetFlags.isPat(:)==1) ,i,ii) = NaN;
                       NetFlags.z(find(NetFlags.z(:,i,ii) < meanP-get(NetFlags.Outlier,'Value')*stdP & NetFlags.isPat(:)==1) ,i,ii) = NaN;
                    end
                end
                cnt = 0;
                %spm('alert"',{[int2str((sum(isnan(NetFlags.z(:)))/numel(NetFlags.z))*100) ' % FC outlieres excluded']});
            else
            %spm('alert"',{'No outliers removed!'});                
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
        if get(NetFlags.Type,'value')>2
            spm('alert*',{'No Monte-Carlo for baseline tests'})
        else
            se_SingleGroup
        end
        
        
    case 'radius'
        NetFlags.Radius = spm_input('VOI Radius','+0','r',NetFlags.Radius,1,[2,10]);
        fg = spm_figure('GetWin','Graphics');
        uicontrol(fg,'Style','Pushbutton','Units','normalized','Position',[0.78 .84 .2 .04],'String',{['Sphere: ' num2str(NetFlags.Radius) 'mm radius']},...
            'HorizontalAlignment','center','FontWeight','bold','Callback','se_Network(''radius'')','ForegroundColor','k');
        
        
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
        weiter =  .06+.72-.07-.025*numel(NetFlags.Label);
        fg = spm_figure('GetWin','Graphics');
        uicontrol(fg,'Style','Pushbutton','Units','normalized','Position',[0.55 weiter-.08 .4 .04],'String',{['p < ' num2str(NetFlags.pvalue)]},...
            'HorizontalAlignment','center','FontWeight','bold','Callback','se_Network(''pvalue'')','ForegroundColor','k');
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
        
        tmp = load(NetFlags.Group);

        if get(NetFlags.PCA,'value')==3
            incl=find(tmp.hasRSfix==1);
        else
            incl=find(tmp.hasRS==1);
        end

        NetFlags.Covariates = tmp.Covariates;
        NetFlags.sub = tmp.sub(incl);
        NetFlags.Cov = tmp.Cov(incl,:);
        NetFlags.isPat = tmp.isPat(incl);
        NetFlags.SubDir = tmp.SubDir(incl);
        try NetFlags.fix = tmp.hasRSfix(incl); end
        try NetFlags.matchQ=tmp.matchQ(incl); end
        
        if sum(ok)==2

            
            GM   = {'NoGMmask','IndGMmask','GroupGMmask'};
            PCA  = {'PCA','NoPCA','FIX'};
            Glob = {'TSR','lTSR','WMCSF','lWMCSF','GSR','lGSR','NoGSR'};
            Filt = {'BP','HP'};
            if NetFlags.Coordinates==0
                des = [GM{get(NetFlags.Mask,'value')} '_' PCA{get(NetFlags.PCA,'value')} '_'...
                    Glob{get(NetFlags.Glob,'value')} '_' Filt{get(NetFlags.Filt,'value')}];
            else
                des = [GM{get(NetFlags.Mask,'value')} '_' PCA{get(NetFlags.PCA,'value')}  '_'...
                    Glob{get(NetFlags.Glob,'value')} '_' Filt{get(NetFlags.Filt,'value')} '_' strrep(num2str(NetFlags.Radius,'%2.1f'),'.','')];
            end
         
            try
                load(fullfile(pwd,'Projects','RS','Networks',spm_str_manip(NetFlags.Group,'rt'),[NetFlags.project '_' des],[NetFlags.project '_' des '.mat']))
                NetFlags.z = z;
                NetFlags.origCorr = z;
                NetFlags.origRanC = zr;
                NetFlags.origParC = pz;
                NetFlags.origRegPC = rpz;

                se_Network('init')
                try; se_Network('adjust'); end 
                
            catch
                se_precomputeNet(des)
                se_Network('init')
                se_Network('adjust')
            end
        end
        
        
        
    case 'mprecomp'
        mgroup = spm_select(Inf,'any',['Select group descriptions'],[],fullfile(pwd,'Lookups'),'.mat',1);
        mvoi   = spm_select(Inf,'any',['Select coordinate list'],[],fullfile(pwd,'VOIs'),'.txt',1);
        
        for gg = 1:size(mgroup,1)
            for vv = 1:size(mvoi,1)
                
                
                NetFlags.Group = deblank(mgroup(gg,:));
                tmp = load(NetFlags.Group);
                
                if get(NetFlags.PCA,'value')==3
                    incl=find(tmp.hasRSfix==1);
                else
                    incl=find(tmp.hasRS==1);
                end
                
                NetFlags.Group = deblank(mgroup(gg,:));
                NetFlags.Covariates = tmp.Covariates;
                NetFlags.sub = tmp.sub(incl);
                NetFlags.Cov = tmp.Cov(incl,:);
                NetFlags.isPat = tmp.isPat(incl);
                NetFlags.SubDir = tmp.SubDir(incl);
                NetFlags.fix = tmp.hasRSfix(incl);
                try NetFlags.matchQ=tmp.matchQ(incl); end
                
                clear tmp
                
                TxtName = deblank(mvoi(vv,:));
                try
                    [NetFlags.X NetFlags.Y NetFlags.Z NetFlags.Label]     = textread(TxtName,'%f %f %f %s');
                    NetFlags.Coordinates = 1;
                catch
                    [NetFlags.VOI NetFlags.Label] = textread(TxtName,'%s %s');
                    NetFlags.Coordinates = 0;
                end
                NetFlags.project = strrep(strrep(spm_str_manip(TxtName,'rt'),'_VOIs',''),'_ImageList','');
                NetFlags.Nodes = zeros(numel(NetFlags.Label,2));
                if ~isfield(NetFlags,'Radius'); NetFlags.Radius = 5; end
                
                for fi =1:2
                    for pc = 1:2
                        for gmm = 1:3
                            for gs = 1:7
                                set(NetFlags.Mask,'value',gmm); set(NetFlags.PCA,'value',pc); set(NetFlags.Glob,'value',gs);
                                GM   = {'NoGMmask','IndGMmask','GroupGMmask'};
                                PCA  = {'PCA','NoPCA','FIX'};
                                Glob = {'TSR','lTSR','WMCSF','lWMCSF','GSR','lGSR','NoGSR'};
                                Filt = {'BP','HP'};

                                if NetFlags.Coordinates==0
                                    des = [GM{get(NetFlags.Mask,'value')} '_' PCA{get(NetFlags.PCA,'value')} '_'...
                                        Glob{get(NetFlags.Glob,'value')} '_' Filt{get(NetFlags.Filt,'value')}];
                                else
                                    des = [GM{get(NetFlags.Mask,'value')} '_' PCA{get(NetFlags.PCA,'value')}  '_'...
                                        Glob{get(NetFlags.Glob,'value')} '_' Filt{get(NetFlags.Filt,'value')} '_'...
                                        strrep(num2str(NetFlags.Radius,'%2.1f'),'.','')];
                                end

                                if exist(fullfile(pwd,'Projects','RS','Networks',spm_str_manip(NetFlags.Group,'rt'),[NetFlags.project '_' des],...
                                        [NetFlags.project '_' des '.mat']))==0
                                    se_precomputeNet(des)
                                end

                            end
                        end
                    end
                end
                         
            end
        end
        
        
    case 'gselect'
        NetFlags.Group = spm_select(1,'any',['Select group description'],[],fullfile(pwd,'Lookups'),'.mat',1);
        tmp = load(NetFlags.Group);
        NetFlags.Covariates = tmp.Covariates;
        
        if get(NetFlags.PCA,'value')==3
            incl=find(tmp.hasRSfix==1);
        else
            incl=find(tmp.hasRS==1);
        end

        NetFlags.sub = tmp.sub(incl);
        NetFlags.Cov = tmp.Cov(incl,:);
        NetFlags.rawCov = tmp.Cov(incl,:);
        NetFlags.isPat = tmp.isPat(incl);
        NetFlags.SubDir = tmp.SubDir(incl);
        try NetFlags.fix = tmp.hasRSfix(incl);
        catch; NetFlags.fix = {}; end
        try NetFlags.matchQ = tmp.matchQ(incl);
        catch; NetFlags.matchQ =[]; end

        clear tmp
        if isfield(NetFlags,'z'); NetFlags = rmfield(NetFlags,'z'); end
        try; NetFlags = rmfield(NetFlags,'activeCov'); end
        try; NetFlags = rmfield(NetFlags,'StatOn');end
        try; NetFlags = rmfield(NetFlags,'Connections');end
        try; NetFlags = rmfield(NetFlags,'Results');end
        try; NetFlags = rmfield(NetFlags,'Description');end
        try; NetFlags = rmfield(NetFlags,'Type');end
        
        try, NetFlags = rmfield(NetFlags,'z'); end
        se_Network('init')
        
        
        
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
        se_Network('init')
        
        
    case 'exit'
        delete(gcf)
        clear
        OverGUI;
        
        
        
    case 'savefiltered'
        
        [file, pfad] = uiputfile({'*.txt','Textfile (*.txt)'},'Save adjusted data as text in',fullfile(pwd,'Projects','VBM','Regional',spm_str_manip(NetFlags.Group,'rt')));
        
        wrz = [];
        wn  = {};
        cnt = 1;
        
        for i=1:numel(NetFlags.Label)-1
            for ii=i+1:numel(NetFlags.Label)
                wrz(cnt,:) = NetFlags.z(:,i,ii)';
                wn{cnt} = [NetFlags.Label{i} ' <-> ' NetFlags.Label{ii}];
                cnt = cnt+1;
            end
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
%         se_figure('Clear','Graphics');
        se_Network('init')
        
end

