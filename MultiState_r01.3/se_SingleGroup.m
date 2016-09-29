function se_SingleGroup

global NetFlags

NetFlags.StatOn = 'S';


try
    Sing = get(NetFlags.Sing,'value');
catch
    Sing = 1;
end

bufferZ = NetFlags.z;

if Sing==2
    NetFlags.z = NetFlags.z(NetFlags.isPat<.5,:,:);
elseif Sing==3
    NetFlags.z = NetFlags.z(NetFlags.isPat>.5,:,:);
end

F = findall(allchild(0),'Flat','Tag','Show Effect');
close(F)
if isfield(NetFlags,'Stats'); NetFlags = rmfield(NetFlags,'Stats'); end
toDo = [];
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


if ~isempty(toDo)
  for i=1:size(toDo,1)
    switch get(NetFlags.Type,'value')
      case 1
        [h rawP(i)] = ttest(NetFlags.z(:,toDo(i,1),toDo(i,2)));
      case 2
        rawP(i)     = signtest(NetFlags.z(:,toDo(i,1),toDo(i,2)));
    end
  end
  
  [rawP I] = sort(rawP,'ascend');
  toDo     = toDo(I,:);

  if Sing==1
      werwar = 'All: ';
  elseif Sing==2
      werwar = 'Controls: ';
  elseif Sing==3
      werwar = 'Patients: ';
  end

    tests = {'T-Tests','Sign-Test'};
  switch get(NetFlags.Method,'value')
    case 1
      pThr    = NetFlags.pvalue;
      desc    = {[werwar  tests{get(NetFlags.Type,'Value')}];['p < ' num2str(NetFlags.pvalue) ' (uncorr.)']};
    case 2
      pThr    = NetFlags.pvalue/size(toDo,1);
      desc    = {[werwar  tests{get(NetFlags.Type,'Value')}];['p < ' num2str(NetFlags.pvalue) ' (Bonf. corr.)']};
    case 3
      pThr    = spm_uc_FDR(NetFlags.pvalue,1,'P',1,sort(rawP,'ascend')',[]);
      desc    = {[werwar  tests{get(NetFlags.Type,'Value')}];['p < ' num2str(NetFlags.pvalue) ' (FDR corr.)']};
  end
  
  toDo = toDo(rawP<pThr,:);
  rawP = rawP(rawP<pThr);
  Results = {};
  if any(toDo)
    for i=1:size(toDo,1)
      Results{i} = [NetFlags.Label{toDo(i,1)} ' <-> ' NetFlags.Label{toDo(i,2)} ' (Z=' num2str(spm_invNcdf(1-rawP(i),0,1),'%3.2f') ')'];
    end
  else
    Results = {'No significant effect'};
  end
  
  NetFlags.Stats.Connections = toDo;
  NetFlags.Stats.Results     = Results;
  NetFlags.Stats.Description = desc;
  NetFlags.Stats.Type        = 'S';
  
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
    'HorizontalAlignment','center','FontSize',FS(12),'FontWeight','bold','BackGroundColor','w');
  
  uicontrol(NetFlags.Stats.fg,'Style','Pushbutton','Units','normalized','Position',[0.55 .955 .4 .04],'String',{'Print all'},...
    'HorizontalAlignment','center','FontWeight','demi','Callback','se_printNet','ForegroundColor','k');
  
  if strcmp(NetFlags.Stats.Results{1}(1:6),'No sig')
    uicontrol(NetFlags.Stats.fg,'Style','text','Units','normalized','Position',[0.1  .85 .8 .04],'String',{NetFlags.Stats.Results{1}},...
      'HorizontalAlignment','center','FontSize',FS(12),'FontAngle','italic','BackGroundColor','w');
  else
    cnt = 1;
    for i=1:numel(NetFlags.Stats.Results)
      toCall = ['se_ShowEffect(' int2str(i) ')'];
      if i<=ceil(numel(NetFlags.Stats.Results)/2)
      uicontrol(NetFlags.Stats.fg,'Style','Pushbutton','Units','normalized','Position',[0.025 .93-(i*.028) .45 .025],'String',{NetFlags.Stats.Results{i}},...
        'HorizontalAlignment','center','FontWeight','demi','Callback',toCall,'ForegroundColor','k','FontSize',FS(8));
      cnt = cnt+1;
      else
      uicontrol(NetFlags.Stats.fg,'Style','Pushbutton','Units','normalized','Position',[0.525 .93-((1+i-cnt)*.028) .45 .025],'String',{NetFlags.Stats.Results{i}},...
        'HorizontalAlignment','center','FontWeight','demi','Callback',toCall,'ForegroundColor','k','FontSize',FS(8));
      end
    end
  end
end

NetFlags.z = bufferZ;
