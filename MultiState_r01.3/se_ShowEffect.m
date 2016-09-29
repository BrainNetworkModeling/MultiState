function se_ShowEffect(choice)

global NetFlags

a = get(0,'ScreenSize');
FS   = spm('FontSizes');
PF   = spm_platform('fonts');
F = findall(allchild(0),'Flat','Tag','Show Effect');
anderes = get(findall(allchild(0),'Flat','Tag','Results'),'Position');
anderes = anderes(1)+anderes(3);
scr = get(0,'ScreenSize');
if isempty(F);
  NetFlags.Stats.se =  figure(...
    'Tag','Show Effect',                          'Position',[min(anderes+20,scr(3)*.6) a(4)*.20 a(4)*.7 a(4)*.9],...
    'Resize','on',                            'MenuBar','figure',                 'Color','w',                              'ColorMap',gray(64),...
    'DefaultTextColor','k',                   'DefaultTextInterpreter','none',     'DefaultTextFontName',PF.helvetica,       'DefaultTextFontSize',FS(12),...
    'DefaultAxesColor','w',                   'DefaultAxesXColor','k',            'DefaultAxesYColor','k',                  'DefaultAxesZColor','k',...
    'DefaultAxesFontName',PF.helvetica,       'DefaultPatchFaceColor','k',        'DefaultPatchEdgeColor','k',              'DefaultSurfaceEdgeColor','k',...
    'DefaultLineColor','k',                   'DefaultUicontrolFontName',PF.helvetica,'NumberTitle','off',...
    'DefaultUicontrolFontSize',FS(12),        'DefaultUicontrolInterruptible','on', 'PaperType','A4',                         'PaperUnits','normalized',...
    'PaperPosition',[.0726 .0644 .854 .870],  'InvertHardcopy','off',             'Renderer','zbuffer',                     'Visible','on','Name','Results details');
else
  set(0,'CurrentFigure',F);
  figure(F); clf
end


switch NetFlags.Stats.Type
  case 'S'
    bar(NetFlags.z(:,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2)))
    set(gca,'XTick',[],'XLim',[0 size(NetFlags.z,1)+1]);
    xlabel('Subjects','FontSize',FS(14))
    ylabel('Functional connectivity (Fishers Z)','FontSize',FS(14))
    title({NetFlags.Stats.Results{choice}; ['in ' spm_str_manip(NetFlags.Group,'rt') ': ' NetFlags.project]},'FontSize',FS(12))
    
  case 'G'
    
      
      if any(NetFlags.isPat)>0
          nowCompare = NetFlags.isPat;
          NetFlags.nowCompareInfo = '';
      else
          nowCompare = NetFlags.nowCompare;
      end

      if size(NetFlags.checkBox,2)>1
          if get(NetFlags.Type,'value') == 1 | get(NetFlags.Type,'value') == 3
              [h p] = ttest(NetFlags.z(nowCompare>0,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2))); if p<0.05; pst = '*'; else, pst = ''; end
              [h p] = ttest(NetFlags.z(nowCompare<1,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2))); if p<0.05; pco = '*'; else, pco = ''; end
              pst = [' (' num2str(nanmean(NetFlags.z(nowCompare>0,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2))),'%4.3f') ')' pst];
              pco = [' (' num2str(nanmean(NetFlags.z(nowCompare<1,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2))),'%4.3f') ')' pco];
          else
              p = signtest(NetFlags.z(nowCompare>0,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2))); if p<0.05; pst = '*'; else, pst = ''; end
              p = signtest(NetFlags.z(nowCompare<1,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2))); if p<0.05; pco = '*'; else, pco = ''; end
              pst = [' (' num2str(median(NetFlags.z(nowCompare>0,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2))),'%4.3f') ')' pst];
              pco = [' (' num2str(median(NetFlags.z(nowCompare<1,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2))),'%4.3f') ')' pco];
          end
      else
      if get(NetFlags.Type,'value') == 1 | get(NetFlags.Type,'value') == 3
          pst = [' (' num2str(nanmean(NetFlags.z(nowCompare>0,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2)))) ')' ];          
          pco = [' (' num2str(nanmean(NetFlags.z(nowCompare<1,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2)))) ')' ];
      else
          pst = [' (' num2str(median(NetFlags.z(nowCompare>0,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2)))) ')' ];          
          pco = [' (' num2str(median(NetFlags.z(nowCompare<1,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2)))) ')' ];
      end
      end
      

    scatter(ones(1,sum(nowCompare>0.5)),(NetFlags.z(nowCompare>0.5,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2))),'MarkerEdgeColor','k','MarkerFaceColor','r','LineWidth',0.5); hold on
    scatter(ones(1,sum(nowCompare<0.5))*2,(NetFlags.z(nowCompare<0.5,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2))),'MarkerEdgeColor','k','MarkerFaceColor','g','LineWidth',0.5); hold on
    if get(NetFlags.Type,'value') == 1  | get(NetFlags.Type,'value') == 3
      scatter([1 2],[nanmean(NetFlags.z(nowCompare>0.5,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2))) nanmean(NetFlags.z(nowCompare<0.5,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2)))],[50 50],'k','filled'); hold on
      line([.9 1.1],[nanmean(NetFlags.z(nowCompare>0.5,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2))) nanmean(NetFlags.z(nowCompare>0.5,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2)))],'Color','k','LIneWidth',1)
      line([1.9 2.1],[nanmean(NetFlags.z(nowCompare<0.5,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2))) nanmean(NetFlags.z(nowCompare<0.5,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2)))],'Color','k','LIneWidth',1)
    else
      scatter([1 2],[nanmedian(NetFlags.z(nowCompare>0.5,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2))) nanmedian(NetFlags.z(nowCompare<0.5,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2)))],[50 50],'k','filled'); hold on
      line([.9 1.1],[nanmedian(NetFlags.z(nowCompare>0.5,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2))) nanmedian(NetFlags.z(nowCompare>0.5,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2)))],'Color','k','LIneWidth',1)
      line([1.9 2.1],[nanmedian(NetFlags.z(nowCompare<0.5,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2))) nanmedian(NetFlags.z(nowCompare<0.5,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2)))],'Color','k','LIneWidth',1)
    end
    if get(NetFlags.Type,'value') == 1 | get(NetFlags.Type,'value') == 3
        d = (nanmean(NetFlags.z(nowCompare>0.5,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2)))-nanmean(NetFlags.z(nowCompare<0.5,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2)))) / sqrt(.5*((nanvar(NetFlags.z(nowCompare>0.5,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2)))+nanvar(NetFlags.z(nowCompare<0.5,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2))))));
      if nanmean(NetFlags.z(nowCompare>0.5,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2))) > nanmean(NetFlags.z(nowCompare<0.5,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2)))
        direction = ['Increase in Patients (d'' = ' num2str(d,'%3.2f') ')'];
      else
        direction = ['Decrease in Patients (d'' = ' num2str(d,'%3.2f') ')'];
      end
    else
        d = (nanmedian(NetFlags.z(nowCompare>0.5,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2)))-nanmedian(NetFlags.z(nowCompare<0.5,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2)))) / sqrt(.5*((nanvar(NetFlags.z(nowCompare>0.5,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2)))+nanvar(NetFlags.z(nowCompare<0.5,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2))))));
      if nanmedian(NetFlags.z(nowCompare>0.5,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2))) > nanmedian(NetFlags.z(nowCompare<0.5,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2)))
        direction = ['Increase in Patients (d'' = ' num2str(d,'%3.2f') ')'];
      else
        direction = ['Decrease in Patients (d'' = ' num2str(d,'%3.2f') ')'];
      end
    end
    if size(NetFlags.checkBox,2)>1;
        ylabel('Functional connectivity (Fishers Z)','FontSize',FS(14));
    else 
        V1 = {'Grey Matter','White Matter'};
        V2 = {'Modulated','Non-linearly modulated','Smoothed Modulated','Smoothed non-linearly modulated'};
        ylabel([V1{get(NetFlags.Mask,'value')} ' ' V2{get(NetFlags.mtype,'value')}],'FontSize',FS(14))
    end
    xlabel(direction,'FontSize',FS(14))
    set(gca,'XLim',[0.5 2.5],'XTick',[1 2],'XTickLabel',{['Patients' pst],['Controls' pco]},'FontSize',FS(12))
    title({NetFlags.Stats.Results{choice}; ['in ' spm_str_manip(NetFlags.Group,'rt') ': ' NetFlags.project]; NetFlags.nowCompareInfo},'FontSize',FS(12))
    
              

    
  case 'C'
    if nargin<1
      choice = get(NetFlags.pop,'value');
    end
    if get(NetFlags.welcheCorr,'value')<3 | get(NetFlags.wasCorr,'Value')>4
      if any(NetFlags.isPat)
        if get(NetFlags.wasCorr,'Value') == 2
        scatter(NetFlags.Cov(NetFlags.isPat>0.5,NetFlags.Stats.Connections(choice,3)),...
          (NetFlags.z(NetFlags.isPat>0.5,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2))),...
          ones(1,sum(NetFlags.isPat>0.5))*40,repmat([.2 .2 .2],sum(NetFlags.isPat>0.5),1),'MarkerEdgeColor','k','MarkerFaceColor','g','LineWidth',0.3); hold on
            xlabel([NetFlags.Covariates{NetFlags.Stats.Connections(choice,3)} '  (only Patients)'],'FontSize',FS(14))
        elseif get(NetFlags.wasCorr,'Value') == 3 | get(NetFlags.wasCorr,'Value') == 1
            scatter(NetFlags.Cov(NetFlags.isPat>0.5,NetFlags.Stats.Connections(choice,3)),...
                (NetFlags.z(NetFlags.isPat>0.5,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2))),...
                ones(1,sum(NetFlags.isPat>0.5))*40,repmat([0 1 0],sum(NetFlags.isPat>0.5),1),'MarkerEdgeColor','k','MarkerFaceColor','g','LineWidth',0.3); hold on
            scatter(NetFlags.Cov(NetFlags.isPat<0.5,NetFlags.Stats.Connections(choice,3)),...
                (NetFlags.z(NetFlags.isPat<0.5,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2))),...
                ones(1,sum(NetFlags.isPat<0.5))*50,repmat([0 0 1],sum(NetFlags.isPat<0.5),1),'MarkerEdgeColor','k','MarkerFaceColor','b','LineWidth',0.3); hold on
            xlabel([NetFlags.Covariates{NetFlags.Stats.Connections(choice,3)} '  (Patients in green)'],'FontSize',FS(14))
        elseif get(NetFlags.wasCorr,'Value') == 5
            scatter(NetFlags.Cov(NetFlags.isPat>0.5,NetFlags.Stats.Connections(choice,3)),...
                (NetFlags.z(NetFlags.isPat>0.5,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2))),...
                ones(1,sum(NetFlags.isPat>0.5))*40,repmat([0 1 0],sum(NetFlags.isPat>0.5),1),'MarkerEdgeColor','k','MarkerFaceColor','g','LineWidth',0.3); hold on
            scatter(NetFlags.Cov(NetFlags.isPat<0.5,NetFlags.Stats.Connections(choice,3)),...
                (NetFlags.z(NetFlags.isPat<0.5,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2))),...
                ones(1,sum(NetFlags.isPat<0.5))*50,repmat([0 0 1],sum(NetFlags.isPat<0.5),1),'MarkerEdgeColor','k','MarkerFaceColor','b','LineWidth',0.3); hold on
            xlabel([NetFlags.Covariates{NetFlags.Stats.Connections(choice,3)} '  (Patients in green)'],'FontSize',FS(14))

        
        elseif get(NetFlags.wasCorr,'Value') == 4
        scatter(NetFlags.Cov(NetFlags.isPat>0.5,NetFlags.Stats.Connections(choice,3)),...
          (NetFlags.z(NetFlags.isPat>0.5,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2))),...
          ones(1,sum(NetFlags.isPat>0.5))*30,repmat([.2 .2 .2],sum(NetFlags.isPat>0.5),1),'MarkerEdgeColor','k','MarkerFaceColor','g','LineWidth',0.3); hold on
            xlabel([NetFlags.Covariates{NetFlags.Stats.Connections(choice,3)} '  (only Patients)'],'FontSize',FS(14))
        end
        
        if get(NetFlags.welcheCorr,'value') < 3  | get(NetFlags.wasCorr,'Value')>4
          if get(NetFlags.wasCorr,'Value') == 1
            STAT = se_regr(NetFlags.Cov(NetFlags.isPat>-Inf,NetFlags.Stats.Connections(choice,3)),(NetFlags.z(NetFlags.isPat>-Inf,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2))));
            plot(STAT.xtmp,STAT.ystar,'LineWidth',2,'LineStyle','-','Color',[0 0 0])
            plot(STAT.xtmp,STAT.cir,'LineWidth',1,'LineStyle','--','Color',[0 0 0])
            
          elseif get(NetFlags.wasCorr,'Value') == 2
            STAT = se_regr(NetFlags.Cov(NetFlags.isPat>0.5,NetFlags.Stats.Connections(choice,3)),(NetFlags.z(NetFlags.isPat>0.5,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2))));
            plot(STAT.xtmp,STAT.ystar,'LineWidth',2,'LineStyle','-','Color','k')
            plot(STAT.xtmp,STAT.cir,'LineWidth',1,'LineStyle','--','Color','k')
            
          elseif get(NetFlags.wasCorr,'Value') == 3
            STAT = se_regr(NetFlags.Cov(NetFlags.isPat>0.5,NetFlags.Stats.Connections(choice,3)),(NetFlags.z(NetFlags.isPat>0.5,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2))));
            plot(STAT.xtmp,STAT.ystar,'LineWidth',2,'LineStyle','-','Color',[0 1 0])
            plot(STAT.xtmp,STAT.cir,'LineWidth',1,'LineStyle','--','Color',[0 1 0])
            
            STAT = se_regr(NetFlags.Cov(NetFlags.isPat<0.5,NetFlags.Stats.Connections(choice,3)),(NetFlags.z(NetFlags.isPat<0.5,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2))));
            plot(STAT.xtmp,STAT.ystar,'LineWidth',2,'LineStyle','-','Color',[0 0 1])
            plot(STAT.xtmp,STAT.cir,'LineWidth',1,'LineStyle','--','Color',[0 0 1])
            
           elseif get(NetFlags.wasCorr,'Value') == 5
            STAT = se_regr(NetFlags.Cov(NetFlags.isPat>0.5,NetFlags.Stats.Connections(choice,3)),(NetFlags.z(NetFlags.isPat>0.5,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2))));
            plot(STAT.xtmp,STAT.ystar,'LineWidth',2,'LineStyle','-','Color',[0 1 0])
            plot(STAT.xtmp,STAT.cir,'LineWidth',1,'LineStyle','--','Color',[0 1 0])
            
            STAT = se_regr(NetFlags.Cov(NetFlags.isPat<0.5,NetFlags.Stats.Connections(choice,3)),(NetFlags.z(NetFlags.isPat<0.5,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2))));
            plot(STAT.xtmp,STAT.ystar,'LineWidth',2,'LineStyle','-','Color',[0 0 1])
            plot(STAT.xtmp,STAT.cir,'LineWidth',1,'LineStyle','--','Color',[0 0 1])
     
            
          elseif get(NetFlags.wasCorr,'Value') == 4
            
            STAT = se_regr(NetFlags.Cov(NetFlags.isPat>0.5,NetFlags.Stats.Connections(choice,3)),(NetFlags.z(NetFlags.isPat>0.5,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2))));
            plot(STAT.xtmp,STAT.ystar,'LineWidth',2,'LineStyle','-','Color',[0 1 0])
            plot(STAT.xtmp,STAT.cir,'LineWidth',1,'LineStyle','--','Color',[0 1 0])
            if get(NetFlags.Type,'value') == 2  | get(NetFlags.Type,'value') == 4
              line([min(NetFlags.Cov(NetFlags.isPat<0.5,NetFlags.Stats.Connections(choice,3))) max(NetFlags.Cov(NetFlags.isPat<0.5,NetFlags.Stats.Connections(choice,3)))],[nanmedian(NetFlags.z(NetFlags.isPat<0.5,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2))) nanmedian(NetFlags.z(NetFlags.isPat<0.5,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2)))],'LineStyle','--','Color','b'); 
              line([min(NetFlags.Cov(NetFlags.isPat>0.5,NetFlags.Stats.Connections(choice,3))) max(NetFlags.Cov(NetFlags.isPat>0.5,NetFlags.Stats.Connections(choice,3)))],[nanmedian(NetFlags.z(NetFlags.isPat>0.5,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2))) nanmedian(NetFlags.z(NetFlags.isPat>0.5,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2)))],'LineStyle','--','Color','g'); 
            else
              line([min(NetFlags.Cov(NetFlags.isPat<0.5,NetFlags.Stats.Connections(choice,3))) max(NetFlags.Cov(NetFlags.isPat<0.5,NetFlags.Stats.Connections(choice,3)))],[nanmean(NetFlags.z(NetFlags.isPat<0.5,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2))) nanmean(NetFlags.z(NetFlags.isPat<0.5,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2)))],'LineStyle','--','Color','b') 
              line([min(NetFlags.Cov(NetFlags.isPat>0.5,NetFlags.Stats.Connections(choice,3))) max(NetFlags.Cov(NetFlags.isPat>0.5,NetFlags.Stats.Connections(choice,3)))],[nanmean(NetFlags.z(NetFlags.isPat>0.5,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2))) nanmean(NetFlags.z(NetFlags.isPat>0.5,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2)))],'LineStyle','--','Color','g') 
            end
            
          end
        end
        
      else
        scatter(NetFlags.Cov(NetFlags.isPat<0.5,NetFlags.Stats.Connections(choice,3)),...
          (NetFlags.z(NetFlags.isPat<0.5,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2))),...
          ones(1,sum(NetFlags.isPat<0.5))*30,repmat([.4 .4 .4],sum(NetFlags.isPat<0.5),1),'o','filled'); hold on
        
        STAT = se_regr(NetFlags.Cov(NetFlags.isPat<0.5,NetFlags.Stats.Connections(choice,3)),(NetFlags.z(NetFlags.isPat<0.5,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2))));
        plot(STAT.xtmp,STAT.ystar,'LineWidth',2,'LineStyle','-','Color',[0 0 1])
        plot(STAT.xtmp,STAT.cir,'LineWidth',1,'LineStyle','--','Color',[0 0 1])
        
        xlabel([NetFlags.Covariates{NetFlags.Stats.Connections(choice,3)} ],'FontSize',FS(14))
      end
      
      axis tight
    if NetFlags.wasison == 1;
        ylabel('Functional connectivity (Fishers Z)','FontSize',FS(14));
    elseif NetFlags.wasison == 2
        V1 = {'Grey Matter','White Matter'};
        V2 = {'Modulated','Non-linearly modulated','Smoothed Modulated','Smoothed non-linearly modulated'};
        ylabel([V1{get(NetFlags.Mask,'value')} ' ' V2{get(NetFlags.mtype,'value')}],'FontSize',FS(14))
    end
%     if get(NetFlags.wasCorr,'Value') == 5
%         title({[spm_str_manip(NetFlags.Group,'rt') ': ' NetFlags.project]; NetFlags.Stats.Results{choice}},'FontSize',FS(12))
%     else 
         title({[spm_str_manip(NetFlags.Group,'rt') ': ' NetFlags.project]; NetFlags.Stats.Results{choice}},'FontSize',FS(12))
%     end
    
    else
      if get(NetFlags.wasCorr,'Value') == 1
        Q = find(NetFlags.isPat>-Inf);
        Q = Q(~isnan(NetFlags.Cov(Q,NetFlags.Stats.Connections(choice,3))));
        labels = unique(NetFlags.Cov(Q,NetFlags.Stats.Connections(choice,3)));
        Xcov     = NetFlags.Cov(Q,NetFlags.Stats.Connections(choice,3));
        [n xout] = hist(Xcov,unique(Xcov));
        n = n(~isnan(xout));
        xout = xout(~isnan(xout));
        ex = find(n<4);
        xQ  = ones(1,numel(Xcov));
        xQ(isnan(Xcov))=0;
        for xii = 1:numel(ex)
          xQ(find(Xcov==xout(ex(xii)))) = 0;
        end
        coupling = (NetFlags.z(Q,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2)));
        coupling = coupling(xQ>0);
        Xcov = Xcov(xQ>0);
        
        labels = unique(Xcov);
        for i=1:numel(labels)
          Label{i} = [NetFlags.Covariates{NetFlags.Stats.Connections(choice,3)} '=' int2str(labels(i))];
        end
        boxplot(coupling,Xcov,'labels',Label)
        
      elseif get(NetFlags.wasCorr,'Value') == 2
        Q = find(NetFlags.isPat>.5);
        Q = Q(~isnan(NetFlags.Cov(Q,NetFlags.Stats.Connections(choice,3))));
        labels = unique(NetFlags.Cov(Q,NetFlags.Stats.Connections(choice,3)));
        Xcov     = NetFlags.Cov(Q,NetFlags.Stats.Connections(choice,3));
        [n xout] = hist(Xcov,unique(Xcov));
        n = n(~isnan(xout));
        xout = xout(~isnan(xout));
        ex = find(n<4);
        xQ  = ones(1,numel(Xcov));
        xQ(isnan(Xcov))=0;
        for xii = 1:numel(ex)
          xQ(find(Xcov==xout(ex(xii)))) = 0;
        end
        coupling = (NetFlags.z(Q,NetFlags.Stats.Connections(choice,1),NetFlags.Stats.Connections(choice,2)));
        coupling = coupling(xQ>0);
        Xcov = Xcov(xQ>0);
        
        labels = unique(Xcov);
        for i=1:numel(labels)
          Label{i} = [NetFlags.Covariates{NetFlags.Stats.Connections(choice,3)} '=' int2str(labels(i))];
        end
        boxplot(coupling,Xcov,'labels',Label)
      end
      
      ylabel('Functional connectivity (Fishers Z)','FontSize',FS(14))
      title({NetFlags.Stats.Results{choice}; ['in ' spm_str_manip(NetFlags.Group,'rt') ': ' NetFlags.project]},'FontSize',FS(12))
    end
    
end
