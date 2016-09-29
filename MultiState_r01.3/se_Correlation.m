function se_Correlation(varargin)


global NetFlags
clc



NetFlags.StatOn = 'C';

if (nargin==0), Action = 'init'; else, Action = varargin{1}; end


switch lower(Action),
    
    case 'init'
        
        try
            WasWar = get(NetFlags.wasCorr,'value');
        catch
            WasWar = 1;
        end
        try
            WelcheWar = get(NetFlags.welcheCorr,'value');
        catch
            WelcheWar = 2;
        end
        try
            hochWar = get(NetFlags.Esize,'value');
        catch
            hochWar = 2;
        end
        
        try
            NetFlags.welcheCorrWar;
        catch
            NetFlags.welcheCorrWar = 1;
        end
        
        fg = spm_figure('GetWin','Graphics');
        weiter =  .06+.72-.07-.025*numel(NetFlags.Label);
            war = get(NetFlags.Method,'Value');
%             NetFlags.Method = uicontrol(fg,'Style','popupmenu','Units','normalized','Position',[0.55 weiter-.15 .4 .04],'Callback','se_Network(''setstat'')',...
%                 'String',str2mat('uncorrected','Bonferroni corrected','FDR corrected'),'Value',war,'ToolTipString','Correction for multiple comparisons');
        
        F = findall(allchild(0),'Flat','Tag','Show Effect');
        close(F)
        if isfield(NetFlags,'Stats'); NetFlags = rmfield(NetFlags,'Stats'); end
        if ~isfield(NetFlags,'activeCov')
            NetFlags.activeCov = zeros(1,numel(NetFlags.Covariates));
        end
        
        spalten = max(ceil(numel(NetFlags.Covariates)/30),2);
        a = get(0,'ScreenSize');
        FS   = spm('FontSizes');
        PF   = spm_platform('fonts');
        F = findall(allchild(0),'Flat','Tag','Results');
        if isempty(F);
            NetFlags.Stats.cov =  figure(...
                'Tag','Results',                          'Position',[a(3)*.01 a(4)*.06 a(3)*spalten*0.17 a(4)*.87],...
                'Resize','on',                            'MenuBar','figure',                 'Color','w',                              'ColorMap',gray(64),...
                'DefaultTextColor','k',                   'DefaultTextInterpreter','tex',     'DefaultTextFontName',PF.helvetica,       'DefaultTextFontSize',FS(12),...
                'DefaultAxesColor','w',                   'DefaultAxesXColor','k',            'DefaultAxesYColor','k',                  'DefaultAxesZColor','k',...
                'DefaultAxesFontName',PF.helvetica,       'DefaultPatchFaceColor','k',        'DefaultPatchEdgeColor','k',              'DefaultSurfaceEdgeColor','k',...
                'DefaultLineColor','k',                   'DefaultUicontrolFontName',PF.helvetica,...
                'DefaultUicontrolFontSize',FS(12),        'DefaultUicontrolInterruptible','on', 'PaperType','A4',                         'PaperUnits','normalized',...
                'PaperPosition',[.0726 .0644 .854 .870],  'InvertHardcopy','off',             'Renderer','zbuffer',                     'Visible','on','Name','Results');
        else
            set(0,'CurrentFigure',F);
            figure(F); clf
            NetFlags.Stats.cov = F;
        end
        
        uicontrol(NetFlags.Stats.cov,'Style','text','Units','normalized','Position',[0  .95 1 .04],'String',{'Available Covariates'},...
            'HorizontalAlignment','center','FontSize',FS(14),'FontWeight','bold','BackGroundColor','w');
        
        
        if ~any(NetFlags.isPat)
            NetFlags.wasCorr = uicontrol(NetFlags.Stats.cov,'Style','popupmenu','Units','normalized','Position',[0.02 .01 .3 .04],'Callback','se_Correlation(''init'')',...
                'String',str2mat('All subjects'),'Value',WasWar,'FontSize',FS(9));
        else
            NetFlags.wasCorr = uicontrol(NetFlags.Stats.cov,'Style','popupmenu','Units','normalized','Position',[0.02 .01 .3 .04],'Callback','se_Correlation(''init'')',...
                'String',str2mat('All subjects','Only Patients','Group Differences','Patients in Main-Effect','Controls/Patients individually'),'Value',WasWar,'FontSize',FS(9));
        end
        
        
        NetFlags.welcheCorr = uicontrol(NetFlags.Stats.cov,'Style','popupmenu','Units','normalized','Position',[0.02 .05 .3 .04],'Callback','se_Correlation(''change'')',...
            'String',str2mat('Linear','Rank','ANOVA'),'Value',WelcheWar,'FontSize',FS(9));
        
        uicontrol(NetFlags.Stats.cov,'Style','Pushbutton','Units','normalized','Position',[0.68 .01 .3 .04],'String',{'Compute'},...
            'HorizontalAlignment','center','ForegroundColor','r','FontWeight','bold','Callback','se_Correlation(''run'')','ForegroundColor','g','FontSize',FS(9));
        
        uicontrol(NetFlags.Stats.cov,'Style','Pushbutton','Units','normalized','Position',[0.68 .05 .3 .04],'String',{'Clear all'},...
            'HorizontalAlignment','center','ForegroundColor','r','FontWeight','bold','Callback','se_Correlation(''clear'')','ForegroundColor','r','FontSize',FS(9));
        
        if get(NetFlags.welcheCorr,'value')<3 & get(NetFlags.wasCorr,'value')<3
            NetFlags.Esize = uicontrol(NetFlags.Stats.cov,'Style','popupmenu','Units','normalized','Position',[0.34 .05 .3 .04],'Callback','se_Correlation(''change'')',...
                'String',str2mat('All','From small effects','From moderate effects','From large effects'),'Value',hochWar,'FontSize',FS(9));
        end
        
        se_Correlation('wascorr')
        
        
    case 'clear'
        
        Q = find(NetFlags.activeCov)
        for i=1:numel(Q)
            se_SetCov(Q(i))
        end
        NetFlags.activeCov = zeros(1,numel(NetFlags.Covariates));
        try; delete(NetFlags.pop); NetFlags = rmfield(NetFlags,'pop'); end
        F = findall(allchild(0),'Flat','Tag','Show Effect');
        close(F)
        
    case 'change'
        FS   = spm('FontSizes');
        try; delete(NetFlags.pop); NetFlags = rmfield(NetFlags,'pop'); end
        F = findall(allchild(0),'Flat','Tag','Show Effect');
        close(F)
        uicontrol(NetFlags.Stats.cov,'Style','text','Units','normalized','Position',[0  .95 1 .04],'String',{'Available Covariates'},...
            'HorizontalAlignment','center','FontSize',FS(14),'FontWeight','bold','BackGroundColor','w');
        if numel(unique([get(NetFlags.welcheCorr,'value') NetFlags.welcheCorrWar]>2.5))>1
            NetFlags.activeCov = zeros(size(NetFlags.activeCov));
            se_Correlation('wascorr')
        end
        NetFlags.welcheCorrWar = get(NetFlags.welcheCorr,'value');
        se_Correlation('init')
        
    case 'wascorr'
        FS   = spm('FontSizes');
        try; delete(NetFlags.pop); NetFlags = rmfield(NetFlags,'pop'); end
        F = findall(allchild(0),'Flat','Tag','Show Effect');
        close(F)
        uicontrol(NetFlags.Stats.cov,'Style','text','Units','normalized','Position',[0  .95 1 .04],'String',{'Available Covariates'},...
            'HorizontalAlignment','center','FontSize',FS(14),'FontWeight','bold','BackGroundColor','w');
        
        spalten = max(ceil(numel(NetFlags.Covariates)/30),2);
        woSpalten = [0:(1/spalten):1];
        beiCon = any(NetFlags.Cov(NetFlags.isPat<.5,:)>0);
        beiPat = any(NetFlags.Cov(NetFlags.isPat>.5,:)>0);
        for i=1:numel(NetFlags.Covariates)
            toCall = ['se_SetCov(' int2str(i) ')'];
            vollSpalten = ceil(i/30);
            
            switch get(NetFlags.wasCorr,'value')
                case 1
                    if beiCon(i)
                        if NetFlags.activeCov(i) == 1
                            uicontrol(NetFlags.Stats.cov,'Style','Pushbutton','Units','normalized','Position',[0.025+woSpalten(vollSpalten) .97-((1+i-((vollSpalten-1)*30))*.028) (.9/spalten) .025],...
                                'String',{NetFlags.Covariates{i}},'HorizontalAlignment','center','FontWeight','bold','Callback',toCall,'ForegroundColor','r','FontSize',FS(10));
                        else
                            uicontrol(NetFlags.Stats.cov,'Style','Pushbutton','Units','normalized','Position',[0.025+woSpalten(vollSpalten) .97-((1+i-((vollSpalten-1)*30))*.028) (.9/spalten) .025],'String',{NetFlags.Covariates{i}},...
                                'HorizontalAlignment','center','FontWeight','demi','Callback',toCall,'ForegroundColor','k','FontSize',FS(8));
                        end
                    else
                        uicontrol(NetFlags.Stats.cov,'Style','Pushbutton','Units','normalized','Position',[0.025+woSpalten(vollSpalten) .97-((1+i-((vollSpalten-1)*30))*.028) (.9/spalten) .025],'String',{NetFlags.Covariates{i}},...
                            'HorizontalAlignment','center','FontWeight','demi','ForegroundColor',[.7 .7 .7],'BackgroundColor',[.9 .9 .9],'FontSize',FS(8));
                    end
                    war = get(NetFlags.welcheCorr,'value');
                    NetFlags.welcheCorr = uicontrol(NetFlags.Stats.cov,'Style','popupmenu','Units','normalized','Position',[0.02 .05 .3 .04],'Callback','se_Correlation(''change'')',...
                        'String',str2mat('Linear','Rank','ANOVA'),'Value',war,'FontSize',FS(9));
                case 2
                    if beiPat(i)
                        if NetFlags.activeCov(i) == 1
                            uicontrol(NetFlags.Stats.cov,'Style','Pushbutton','Units','normalized','Position',[0.025+woSpalten(vollSpalten) .97-((1+i-((vollSpalten-1)*30))*.028) (.9/spalten) .025],...
                                'String',{NetFlags.Covariates{i}},'HorizontalAlignment','center','FontWeight','bold','Callback',toCall,'ForegroundColor','r','FontSize',FS(10));
                        else
                            uicontrol(NetFlags.Stats.cov,'Style','Pushbutton','Units','normalized','Position',[0.025+woSpalten(vollSpalten) .97-((1+i-((vollSpalten-1)*30))*.028) (.9/spalten) .025],'String',{NetFlags.Covariates{i}},...
                                'HorizontalAlignment','center','FontWeight','demi','Callback',toCall,'ForegroundColor','k','FontSize',FS(8));
                        end
                    else
                        uicontrol(NetFlags.Stats.cov,'Style','Pushbutton','Units','normalized','Position',[0.025+woSpalten(vollSpalten) .97-((1+i-((vollSpalten-1)*30))*.028) (.9/spalten) .025],'String',{NetFlags.Covariates{i}},...
                            'HorizontalAlignment','center','FontWeight','demi','ForegroundColor',[.7 .7 .7],'BackgroundColor',[.9 .9 .9],'FontSize',FS(8));
                    end
                    war = get(NetFlags.welcheCorr,'value');
                    NetFlags.welcheCorr = uicontrol(NetFlags.Stats.cov,'Style','popupmenu','Units','normalized','Position',[0.02 .05 .3 .04],'Callback','se_Correlation(''change'')',...
                        'String',str2mat('Linear','Rank','ANOVA'),'Value',war,'FontSize',FS(9));
                case 3
                    if beiPat(i) & beiCon(i)
                        if NetFlags.activeCov(i) == 1
                            uicontrol(NetFlags.Stats.cov,'Style','Pushbutton','Units','normalized','Position',[0.025+woSpalten(vollSpalten) .97-((1+i-((vollSpalten-1)*30))*.028) (.9/spalten) .025],...
                                'String',{NetFlags.Covariates{i}},'HorizontalAlignment','center','FontWeight','bold','Callback',toCall,'ForegroundColor','r','FontSize',FS(10));
                        else
                            uicontrol(NetFlags.Stats.cov,'Style','Pushbutton','Units','normalized','Position',[0.025+woSpalten(vollSpalten) .97-((1+i-((vollSpalten-1)*30))*.028) (.9/spalten) .025],'String',{NetFlags.Covariates{i}},...
                                'HorizontalAlignment','center','FontWeight','demi','Callback',toCall,'ForegroundColor','k','FontSize',FS(8));
                        end
                    else
                        uicontrol(NetFlags.Stats.cov,'Style','Pushbutton','Units','normalized','Position',[0.025+woSpalten(vollSpalten) .97-((1+i-((vollSpalten-1)*30))*.028) (.9/spalten) .025],'String',{NetFlags.Covariates{i}},...
                            'HorizontalAlignment','center','FontWeight','demi','ForegroundColor',[.7 .7 .7],'BackgroundColor',[.9 .9 .9],'FontSize',FS(8));
                    end
                    if get(NetFlags.welcheCorr,'value')==3
                        NetFlags.welcheCorr = uicontrol(NetFlags.Stats.cov,'Style','popupmenu','Units','normalized','Position',[0.02 .05 .3 .04],'Callback','se_Correlation(''change'')',...
                            'String',str2mat('Linear','Rank'),'Value',2);
                    else
                        war = get(NetFlags.welcheCorr,'value');
                        NetFlags.welcheCorr = uicontrol(NetFlags.Stats.cov,'Style','popupmenu','Units','normalized','Position',[0.02 .05 .3 .04],'Callback','se_Correlation(''change'')',...
                            'String',str2mat('Linear','Rank'),'Value',war,'FontSize',FS(9));
                    end
                case 4
                    if beiPat(i)
                        if NetFlags.activeCov(i) == 1
                            uicontrol(NetFlags.Stats.cov,'Style','Pushbutton','Units','normalized','Position',[0.025+woSpalten(vollSpalten) .97-((1+i-((vollSpalten-1)*30))*.028) (.9/spalten) .025],...
                                'String',{NetFlags.Covariates{i}},'HorizontalAlignment','center','FontWeight','bold','Callback',toCall,'ForegroundColor','r','FontSize',FS(10));
                        else
                            uicontrol(NetFlags.Stats.cov,'Style','Pushbutton','Units','normalized','Position',[0.025+woSpalten(vollSpalten) .97-((1+i-((vollSpalten-1)*30))*.028) (.9/spalten) .025],'String',{NetFlags.Covariates{i}},...
                                'HorizontalAlignment','center','FontWeight','demi','Callback',toCall,'ForegroundColor','k','FontSize',FS(8));
                        end
                    else
                        uicontrol(NetFlags.Stats.cov,'Style','Pushbutton','Units','normalized','Position',[0.025+woSpalten(vollSpalten) .97-((1+i-((vollSpalten-1)*30))*.028) (.9/spalten) .025],'String',{NetFlags.Covariates{i}},...
                            'HorizontalAlignment','center','FontWeight','demi','ForegroundColor',[.7 .7 .7],'BackgroundColor',[.9 .9 .9],'FontSize',FS(8));
                    end
                    war = get(NetFlags.welcheCorr,'value');
                    NetFlags.welcheCorr = uicontrol(NetFlags.Stats.cov,'Style','popupmenu','Units','normalized','Position',[0.02 .05 .3 .04],'Callback','se_Correlation(''change'')',...
                        'String',str2mat('Linear','Rank','ANOVA'),'Value',war,'FontSize',FS(9));
                 case 5
                    if beiPat(i) & beiCon(i)
                        if NetFlags.activeCov(i) == 1
                            uicontrol(NetFlags.Stats.cov,'Style','Pushbutton','Units','normalized','Position',[0.025+woSpalten(vollSpalten) .97-((1+i-((vollSpalten-1)*30))*.028) (.9/spalten) .025],...
                                'String',{NetFlags.Covariates{i}},'HorizontalAlignment','center','FontWeight','bold','Callback',toCall,'ForegroundColor','r','FontSize',FS(10));
                        else
                            uicontrol(NetFlags.Stats.cov,'Style','Pushbutton','Units','normalized','Position',[0.025+woSpalten(vollSpalten) .97-((1+i-((vollSpalten-1)*30))*.028) (.9/spalten) .025],'String',{NetFlags.Covariates{i}},...
                                'HorizontalAlignment','center','FontWeight','demi','Callback',toCall,'ForegroundColor','k','FontSize',FS(8));
                        end
                    else
                        uicontrol(NetFlags.Stats.cov,'Style','Pushbutton','Units','normalized','Position',[0.025+woSpalten(vollSpalten) .97-((1+i-((vollSpalten-1)*30))*.028) (.9/spalten) .025],'String',{NetFlags.Covariates{i}},...
                            'HorizontalAlignment','center','FontWeight','demi','ForegroundColor',[.7 .7 .7],'BackgroundColor',[.9 .9 .9],'FontSize',FS(8));
                    end
                    war = get(NetFlags.welcheCorr,'value');
                    NetFlags.welcheCorr = uicontrol(NetFlags.Stats.cov,'Style','popupmenu','Units','normalized','Position',[0.02 .05 .3 .04],'Callback','se_Correlation(''change'')',...
                        'String',str2mat('Linear','Rank','ANOVA'),'Value',war,'FontSize',FS(9));

            end
        end
        
    case 'run'
        if size(NetFlags.checkBox,2)==2
            toDo = [];
            nowCov = find(NetFlags.activeCov);
            for c = 1:numel(nowCov)
                for i=1:numel(NetFlags.Label)
                    for ii=1:numel(NetFlags.Label)
                        if i~=ii
                            if get(NetFlags.checkBox(i,1),'Value') & get(NetFlags.checkBox(ii,2),'Value')
                                if isempty(toDo) || ~any(sum(toDo == repmat([ii i nowCov(c)],size(toDo,1),1),2)==3)
                                    if get(NetFlags.wasCorr,'Value')<4 | get(NetFlags.wasCorr,'Value')>4
                                        toDo = [toDo; [i ii nowCov(c)]];
                                    else
                                        
                                        switch get(NetFlags.Type,'value')
                                            case 1
                                                [h nowP] = ttest2(NetFlags.z(NetFlags.isPat>.5,i,ii),NetFlags.z(NetFlags.isPat<.5,i,ii));
                                            case 2
                                                nowP     = ranksum(NetFlags.z(NetFlags.isPat>.5,i,ii),NetFlags.z(NetFlags.isPat<.5,i,ii));
                                            case 3
                                                nowP     = se_permute(NetFlags.z(:,i),NetFlags.isPat,0);
                                            case 4
                                                nowP     = se_permute(NetFlags.z(:,i),NetFlags.isPat,1);
                                        end
                                        if nowP < NetFlags.pvalue
                                            toDo = [toDo; [i ii nowCov(c)]];
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
            
        else
            
            toDo = [];
            nowCov = find(NetFlags.activeCov);
            for c = 1:numel(nowCov)
                for i=1:numel(NetFlags.Label)
                    if get(NetFlags.checkBox(i,1),'Value')
                        if get(NetFlags.wasCorr,'Value')<4 | get(NetFlags.wasCorr,'Value')>4
                            toDo = [toDo; [i 1 nowCov(c)]];
                        else
                            switch get(NetFlags.Type,'value')
                                case 1
                                    [h nowP] = ttest2(NetFlags.z(NetFlags.isPat>.5,i),NetFlags.z(NetFlags.isPat<.5,i));
                                case 2
                                    nowP     = ranksum(NetFlags.z(NetFlags.isPat>.5,i),NetFlags.z(NetFlags.isPat<.5,i));
                                case 3
                                    nowP     = se_permute(NetFlags.z(:,i),NetFlags.isPat,0);
                                case 4
                                    nowP     = se_permute(NetFlags.z(:,i),NetFlags.isPat,1);
                            end
                            if nowP < NetFlags.pvalue
                                toDo = [toDo; [i 1 nowCov(c)]];
                            end
                        end
                    end
                end
            end
            
            
            
        end
        
        
        
        
        if ~isempty(toDo)
            
            switch (get(NetFlags.wasCorr,'Value'));
                case 1
                    Q = NetFlags.isPat>-1; singleCor = true;
                case 2
                    Q = NetFlags.isPat>.5; singleCor = true;
                case 3
                    singleCor = false;
                case 4
                    Q = NetFlags.isPat>.5; singleCor = true;
                case 5
                    Q = NetFlags.isPat>.5; singleCor = true;
                    Q2= NetFlags.isPat<.5; singleCor = true;
            end
            
            switch get(NetFlags.welcheCorr,'value')
                case 1
                    corrtype = 'Pearson'; corrF = true;
                case 2
                    corrtype = 'Spearman'; corrF = true;
                case 3
                    corrF = false;
                case 4
                    corrF = true;
            end
            
            if singleCor && corrF
                for i=1:size(toDo,1)
                    [rho(i),rawP(i)] = corr(NetFlags.z(Q,toDo(i,1),toDo(i,2)),NetFlags.Cov(Q,toDo(i,3)),'type',corrtype,'rows','pairwise');
                    if exist('Q2','var')
                        [rho2(i),rawP2(i)] = corr(NetFlags.z(Q2,toDo(i,1),toDo(i,2)),NetFlags.Cov(Q2,toDo(i,3)),'type',corrtype,'rows','pairwise');
                    end 
                end
            elseif ~singleCor  && corrF
                for i=1:size(toDo,1)
                    
                    [r1 p1] = corr(NetFlags.z(NetFlags.isPat>.5,toDo(i,1),toDo(i,2)),NetFlags.Cov(NetFlags.isPat>.5,toDo(i,3)),'type',corrtype,'rows','pairwise');
                    [r2 p2] = corr(NetFlags.z(NetFlags.isPat<.5,toDo(i,1),toDo(i,2)),NetFlags.Cov(NetFlags.isPat<.5,toDo(i,3)),'type',corrtype,'rows','pairwise');
                    
                    Z = (atanh(r1)-atanh(r2))/sqrt((1/(sum(NetFlags.isPat)-3))+(1/(sum(~(NetFlags.isPat))-3)));
                    rawP(i) = min(spm_Ncdf(Z), spm_Ncdf(-Z));
                    if max(p1,p2)>(NetFlags.pvalue*2)
                        rawP(i) = 1;
                    end
                    rho(i) = min(spm_Ncdf(Z));
                end
            elseif ~corrF && (get(NetFlags.wasCorr,'Value'))<3 | get(NetFlags.wasCorr,'Value')>4
                for i=1:size(toDo,1)
                    Xcov     = NetFlags.Cov(Q,toDo(i,3));
                    rho(i) = 0;
                    rawP(i)= 1;
                    if numel(unique(Xcov))<(numel(Xcov)/3)
                        [n xout] = hist(Xcov,unique(Xcov));
                        n = n(~isnan(xout));
                        xout = xout(~isnan(xout));
                        if sum(n>3)>1
                            ex = find(n<4);
                            xQ  = ones(1,numel(Xcov));
                            xQ(isnan(Xcov))=0;
                            
                            for xii = 1:numel(ex)
                                xQ(find(Xcov==xout(ex(xii)))) = 0;
                            end
                            coupling = NetFlags.z(Q,toDo(i,1),toDo(i,2));
                            
                            coupling = coupling(xQ>0);
                            Xcov = Xcov(xQ>0);
                            
                            [n xout] = hist(Xcov,unique(Xcov));
                            [buffer,table,stats] = anova1(coupling,Xcov,'off');
                            if numel(unique(sign(diff(stats.means))))==1
                                rawP(i) = buffer;
                                rho(i) = table{2,end-1};
                            end
                        end
                    end
                end
            end
            
            switch get(NetFlags.Method,'value')
                case 1
                    pThr    = NetFlags.pvalue;
                    desc    = ['p < ' num2str(NetFlags.pvalue) ' (uncorr.)'];
                case 2
                    pThr    = NetFlags.pvalue/size(toDo,1);
                    desc    = ['p < ' num2str(NetFlags.pvalue) ' (Bonf. corr.)'];
                case 3
                    pThr    = max(spm_uc_FDR(NetFlags.pvalue,1,'P',1,sort(rawP,'ascend')',[]),NetFlags.pvalue/size(toDo,1));
                    desc    = ['p < ' num2str(NetFlags.pvalue) ' (FDR corr.)'];
            end
            
            [B I] = sort(rawP,'ascend');
            rawP = rawP(I);
            toDo = toDo(I,:);
            rho  = rho(I);
            
            if exist('Q2','var')
                if sum(rawP<pThr) > sum(rawP2<pThr)
                    rawP2 = rawP2(I);
                    rho2  = rho2(I);                    
                    rawP2 = rawP2(rawP<pThr);
                    rawP = rawP(rawP<pThr);
                    toDo = toDo(rawP<pThr,:);
                    rho2 = rho2(rawP<pThr);
                    rho  = rho(rawP<pThr);    
                else
                    [B I] = sort(rawP2,'ascend');
                    rawP = rawP(I);
                    toDo = toDo(I,:);
                    rho  = rho(I);
                    rawP2 = rawP2(I);
                    rho2  = rho2(I);                    
                    rawP2 = rawP2(rawP2<pThr);
                    rawP = rawP(rawP2<pThr);
                    toDo = toDo(rawP2<pThr,:);
                    rho2 = rho2(rawP2<pThr);
                    rho  = rho(rawP2<pThr);                   
                end    
            else
                rawP = rawP(rawP<pThr);
                toDo = toDo(rawP<pThr,:);
                rho = rho(rawP<pThr);
            end
            
            if singleCor && corrF
                try
                    switch get(NetFlags.Esize,'value');
                        case 2; I = find(abs(rho)>=0.1);
                        case 3; I = find(abs(rho)>=0.3);
                        case 4; I = find(abs(rho)>=0.5);
                        otherwise; I = 1:numel(rho);
                    end
                catch
%                    I = find(abs(rho)>=0.1);
                    I = 1:numel(rho);
                end
                rawP = rawP(I);
                toDo = toDo(I,:);
                rho  = rho(I);
            end
            
            Results = {};
            if any(toDo)
                for i=1:size(toDo,1)
                    if size(NetFlags.checkBox,2)==2
                    if singleCor && corrF && exist('Q2','var')
                        Results{i} = [NetFlags.Label{toDo(i,1)} ' <-> ' NetFlags.Label{toDo(i,2)}...
                            ' ~ ' NetFlags.Covariates{toDo(i,3)},...
                            ' (CON: r=' num2str(rho2(i),'%3.2f') ' / p=' num2str(rawP2(i),'%5.4f')...
                            ' PAT: r=' num2str(rho(i),'%3.2f') ' / p=' num2str(rawP(i),'%5.4f') ')'];
                    elseif singleCor && corrF   
                        Results{i} = [NetFlags.Label{toDo(i,1)} ' <-> ' NetFlags.Label{toDo(i,2)}...
                            ' ~ ' NetFlags.Covariates{toDo(i,3)},...
                            '   (r=' num2str(rho(i),'%3.2f') ' /'...
                            ' p=' num2str(rawP(i),'%5.4f') ')'];
                    elseif ~singleCor  && corrF
                        Results{i} = [NetFlags.Label{toDo(i,1)} ' <-> ' NetFlags.Label{toDo(i,2)}...
                            ' ~ ' NetFlags.Covariates{toDo(i,3)},...
                            '   p=' num2str(rawP(i),'%5.4f')];
                    elseif ~corrF
                        Results{i} = [NetFlags.Label{toDo(i,1)} ' <-> ' NetFlags.Label{toDo(i,2)}...
                            '   ' NetFlags.Covariates{toDo(i,3)},...
                            ': F=' num2str(rho(i),'%3.3f'),...
                            ' (p=' num2str(rawP(i),'%3.3f') ')'];
                    else
                          Results{i} = [NetFlags.Label{toDo(i,1)} ' <-> ' NetFlags.Label{toDo(i,2)}...
                            ' ~ ' NetFlags.Covariates{toDo(i,3)},...
                            ' (CON: r=' num2str(rho(i),'%3.2f') '/' ' p=' num2str(rawP(i),'%5.4f')...
                            ' PAT: r=' num2str(rho2(i),'%3.2f') '/' ' p=' num2str(rawP2(i),'%5.4f') ')'];                      
                    end
                   else
                    if singleCor && corrF && exist('Q2','var')
                        Results{i} = [NetFlags.Label{toDo(i,1)}...
                            ' ~ ' NetFlags.Covariates{toDo(i,3)},...
                            ' (CON: r=' num2str(rho2(i),'%3.2f') ' / p=' num2str(rawP2(i),'%5.4f')...
                            ' PAT: r=' num2str(rho(i),'%3.2f') ' / p=' num2str(rawP(i),'%5.4f') ')'];
                    elseif singleCor && corrF   
                        Results{i} = [NetFlags.Label{toDo(i,1)} ...
                            ' ~ ' NetFlags.Covariates{toDo(i,3)},...
                            '   (r=' num2str(rho(i),'%3.2f') ' /'...
                            '   p=' num2str(rawP(i),'%5.4f') ')'];
                    elseif ~singleCor  && corrF
                        Results{i} = [NetFlags.Label{toDo(i,1)}...
                            ' ~ ' NetFlags.Covariates{toDo(i,3)},...
                            '   p=' num2str(rawP(i),'%5.4f')];
                    elseif ~corrF
                        Results{i} = [NetFlags.Label{toDo(i,1)} ...
                            '   ' NetFlags.Covariates{toDo(i,3)},...
                            ': F=' num2str(rho(i),'%3.3f'),...
                            ' (p=' num2str(rawP(i),'%3.3f') ')'];
                    end
                    end
                end
                
                NetFlags.Stats.Connections = toDo;
                NetFlags.Stats.Results     = Results;
                NetFlags.Stats.Description = desc;
                NetFlags.Stats.Type        = 'C';
                
                
                % Drop-Down
                NetFlags.pop = uicontrol(NetFlags.Stats.cov,'Style','popupmenu','Units','normalized','Position',[0.01  .95 .78 .04],'Callback','se_ShowEffect',...
                    'String',str2mat(Results),'Value',1);
                
                uicontrol(NetFlags.Stats.cov,'Style','Pushbutton','Units','normalized','Position',[0.81 .955 .18 .035],'String',{'Print all'},...
                    'HorizontalAlignment','center','FontWeight','demi','Callback','se_printNet','ForegroundColor','k');
                
            else
                Results = {'No significant effect'};
                NetFlags.pop = uicontrol(NetFlags.Stats.cov,'Style','popupmenu','Units','normalized','Position',[0.01  .95 .78 .04],...
                    'String',str2mat(Results),'Value',1);
            end
            
        end
end


end

