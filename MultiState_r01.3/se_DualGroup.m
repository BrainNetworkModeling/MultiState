function se_DualGroup

global NetFlags

NetFlags.StatOn = 'G';

F = findall(allchild(0),'Flat','Tag','Show Effect');
close(F)
if isfield(NetFlags,'Stats'); NetFlags = rmfield(NetFlags,'Stats'); end

if any(NetFlags.isPat)>0
    nowCompare = NetFlags.isPat;
else
    Q = zeros(1,size(NetFlags.Cov,2));
    for i=1:size(NetFlags.Cov,2)
        if all(isnan(NetFlags.rawCov(:,i))==0) & numel(unique(NetFlags.rawCov(:,i)))==2
            Q(i) = 1;
        end
    end
    if numel(find(Q))>3 | numel(find(Q))<1
        spm('alert*','Could not find anything that may be a group variable'); return
    else
        Q = find(Q); str = [];
        for i=1:numel(Q);
            str = [str strtok(NetFlags.Covariates{Q(i)},' ') '|'];
        end
        str = str(1:end-1);
        what = spm_input('Group by','+1','bd',str,[1:numel(Q)],1);
        groups = unique(NetFlags.Cov(:,Q(what)));
        nowCompare = NetFlags.isPat;
        nowCompare(NetFlags.rawCov(:,Q(what))==groups(1)) = 0;
        nowCompare(NetFlags.rawCov(:,Q(what))==groups(2)) = 1;
        spm('alert!',['Subjects with value  ' int2str(groups(2)) '  for  "' NetFlags.Covariates{Q(what)} '"  are labeled as patients']);
        
        NetFlags.nowCompare = nowCompare;
        NetFlags.nowCompareInfo = ['Subjects with   "' int2str(groups(2)) '"   in   "' NetFlags.Covariates{Q(what)} '"   labeled Patients'];
    end
end

toDo = [];


if size(NetFlags.checkBox,2)>1
    for i=1:numel(NetFlags.Label)
        for ii=1:numel(NetFlags.Label)
            if i~=ii
                if get(NetFlags.checkBox(i,1),'Value') & get(NetFlags.checkBox(ii,2),'Value')
                    if isempty(toDo) || ~any(sum(toDo == repmat([ii i],size(toDo,1),1),2)==2)
                        toDo = [toDo; [i ii]];
                    end
                end
            end
        end
    end
else
    for i=1:numel(NetFlags.Label)
        if get(NetFlags.checkBox(i,1),'Value')
            toDo = [toDo; i 1];
        end
    end
end




if numel(toDo)>0
  
    
  for i=1:size(toDo,1)
    switch get(NetFlags.Type,'value')
      case 1
        [h rawP(i)] = ttest2(NetFlags.z(nowCompare>.5,toDo(i,1),toDo(i,2)),NetFlags.z(nowCompare<.5,toDo(i,1),toDo(i,2)));
      case 2
        rawP(i)     = ranksum(NetFlags.z(nowCompare>.5,toDo(i,1),toDo(i,2)),NetFlags.z(nowCompare<.5,toDo(i,1),toDo(i,2)));
      case 3
        rawP(i)     = se_permute(NetFlags.z(:,toDo(i,1),toDo(i,2)),nowCompare,0);
      case 4
        rawP(i)     = se_permute(NetFlags.z(:,toDo(i,1),toDo(i,2)),nowCompare,1);
    end
  end
  
  
  
  
  [rawP I] = sort(rawP,'ascend');
  toDo     = toDo(I,:);
  
  tests = {'T-Tests','Ranksum','MC (Mean)','MC (Median)'};
  switch get(NetFlags.Method,'value')
      case 1
          pThr    = NetFlags.pvalue;
          desc    = {['Group comparison by '  tests{get(NetFlags.Type,'Value')}];['p < ' num2str(NetFlags.pvalue) ' (uncorr.)']};
      case 2
          pThr    = NetFlags.pvalue/size(toDo,1);
          desc    = {['Group comparison by '  tests{get(NetFlags.Type,'Value')}];['p < ' num2str(NetFlags.pvalue) ' (Bonf. corr.)']};
      case 3
          pThr    = max(spm_uc_FDR(NetFlags.pvalue,1,'P',1,sort(rawP,'ascend')',[]),NetFlags.pvalue/size(toDo,1));
          desc    = {['Group comparison by '  tests{get(NetFlags.Type,'Value')}];['p < ' num2str(NetFlags.pvalue) ' (FDR corr.)']};
  end
  
  
  toDo = toDo(rawP<pThr,:);
  Results = {};
  if any(toDo)
      for i=1:size(toDo,1)
          if get(NetFlags.Type,'value')<3
              if size(NetFlags.checkBox,2)>1
                  Results{i} = [NetFlags.Label{toDo(i,1)} ' <-> ' NetFlags.Label{toDo(i,2)} ' (p<' num2str(max(rawP(i),0.0001),'%5.4f') ')'];
              else
                  Results{i} = [NetFlags.Label{toDo(i,1)} '(p=' num2str(max(rawP(i),0.0001),'%5.4f') ')'];
              end
          else
              if size(NetFlags.checkBox,2)>1
                  Results{i} = [NetFlags.Label{toDo(i,1)} ' <-> ' NetFlags.Label{toDo(i,2)} ' (P=' num2str((1-max(rawP(i),0.0001)),'%3.2f') ')'];
              else
                  Results{i} = [NetFlags.Label{toDo(i,1)} ' (P=' num2str((1-max(rawP(i),0.0001)),'%3.2f') ')'];
              end
          end
      end
  else
    Results = {'No significant effect'};
  end
  
  NetFlags.Stats.Connections = toDo;
  NetFlags.Stats.Results     = Results;
  NetFlags.Stats.Description = desc;
  NetFlags.Stats.Type        = 'G';
  
  a = get(0,'ScreenSize');
  FS   = spm('FontSizes');
  PF   = spm_platform('fonts');
  F = findall(allchild(0),'Flat','Tag','Results');
  if isempty(F);
    NetFlags.Stats.fg =  figure(...
      'Tag','Results',                          'Position',[a(3)*.01 a(4)*.05 a(3)*.3 a(4)*.9],...
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
    NetFlags.Stats.fg = F;
  end
  
  uicontrol(NetFlags.Stats.fg,'Style','text','Units','normalized','Position',[0  .94 .5 .05],'String',NetFlags.Stats.Description,...
    'HorizontalAlignment','center','FontSize',12,'FontWeight','bold','BackGroundColor','w');

  uicontrol(NetFlags.Stats.fg,'Style','Pushbutton','Units','normalized','Position',[0.55 .955 .4 .04],'String',{'Print all'},...
    'HorizontalAlignment','center','FontWeight','demi','Callback','se_printNet','ForegroundColor','k');
  
  if strcmp(NetFlags.Stats.Results{1}(1:3),'No ')
    uicontrol(NetFlags.Stats.fg,'Style','text','Units','normalized','Position',[0.1  .85 .8 .04],'String',{NetFlags.Stats.Results{1}},...
      'HorizontalAlignment','center','FontSize',FS(12),'FontAngle','italic','BackGroundColor','w');
  else
    cnt = 1;
    for i=1:numel(NetFlags.Stats.Results)
      toCall = ['se_ShowEffect(' int2str(i) ')'];
      if i<=ceil(numel(NetFlags.Stats.Results)/2)
      uicontrol(NetFlags.Stats.fg,'Style','Pushbutton','Units','normalized','Position',[0.015 .93-(i*.028) .475 .025],'String',{NetFlags.Stats.Results{i}},...
        'HorizontalAlignment','center','FontWeight','demi','Callback',toCall,'ForegroundColor','k','FontSize',FS(8));
      cnt = cnt+1;
      else
      uicontrol(NetFlags.Stats.fg,'Style','Pushbutton','Units','normalized','Position',[0.515 .93-((1+i-cnt)*.028) .475 .025],'String',{NetFlags.Stats.Results{i}},...
        'HorizontalAlignment','center','FontWeight','demi','Callback',toCall,'ForegroundColor','k','FontSize',FS(8));
      end
    end
  end
end