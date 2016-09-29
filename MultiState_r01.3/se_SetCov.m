function se_SetCov(i)

global NetFlags

spalten = max(ceil(numel(NetFlags.Covariates)/30),2);
woSpalten = [0:(1/spalten):1];
toCall = ['se_SetCov(' int2str(i) ')'];
vollSpalten = ceil(i/30);
FS   = spm('FontSizes');

try; delete(NetFlags.pop); NetFlags = rmfield(NetFlags,'pop'); end
F = findall(allchild(0),'Flat','Tag','Show Effect');
close(F)
uicontrol(NetFlags.Stats.cov,'Style','text','Units','normalized','Position',[0.05  .95 .9 .04],'String',{'Available Covariates'},...
  'HorizontalAlignment','center','FontSize',FS(14),'FontWeight','bold','BackGroundColor','w');

if get(NetFlags.welcheCorr,'value')<3
  switch(get(NetFlags.wasCorr,'Value'));
    case 1
      Q = find(NetFlags.isPat>-1);
    case 2
      Q = find(NetFlags.isPat>.5);
    case 3
      Q = find(NetFlags.isPat>-Inf);
    case 4
      Q = find(NetFlags.isPat>.5);
  end
  if 0 % NetFlags.activeCov(i) == 0 && all(round(NetFlags.Cov(Q(~isnan(NetFlags.Cov(Q,i))),i))==(NetFlags.Cov(Q(~isnan(NetFlags.Cov(Q,i))),i))) && numel(unique(NetFlags.Cov(Q(~isnan(NetFlags.Cov(Q,i))),i)))<4
    sure = questdlg('Sure this is NOT categorical data?','CAVE!','Yes','No','No');
    if sure(1) == 'Y';
      OK = 1;
    else
      OK = 0;
    end
  else
    OK = 1;
  end
end


if get(NetFlags.welcheCorr,'value')==3
  switch(get(NetFlags.wasCorr,'Value'));
    case 1
      Q = find(NetFlags.isPat>-1);
    case 2
      Q = find(NetFlags.isPat>.5);
  end
  if 0 % NetFlags.activeCov(i) == 0 && any(round(NetFlags.Cov(Q(~isnan(NetFlags.Cov(Q,i))),i))~=(NetFlags.Cov(Q(~isnan(NetFlags.Cov(Q,i))),i))) || numel(unique(NetFlags.Cov(Q(~isnan(NetFlags.Cov(Q,i))),i)))>5
    sure = questdlg('Sure this is categorical data?','CAVE!','Yes','No','No');
    if sure(1) == 'Y';
      OK = 1;
    else
      OK = 0;
    end
  else
    OK = 1;
  end
end


if OK
  if NetFlags.activeCov(i) == 0
    uicontrol(NetFlags.Stats.cov,'Style','Pushbutton','Units','normalized','Position',[0.025+woSpalten(vollSpalten) .97-((1+i-((vollSpalten-1)*30))*.028) (.9/spalten) .025],...
      'String',{NetFlags.Covariates{i}},'HorizontalAlignment','center','FontWeight','bold','Callback',toCall,'ForegroundColor','r','FontSize',FS(10));
    NetFlags.activeCov(i) = 1;
    
%     if get(NetFlags.welcheCorr,'value')==3
%       if sum(NetFlags.activeCov)>1
%         Q = NetFlags.activeCov;
%         Q(i) = 0;
%         Q = find(Q);
%         for xi = 1:numel(Q)
%           i = Q(xi);
%           toCall = ['se_SetCov(' int2str(i) ')'];
%           
%           uicontrol(NetFlags.Stats.cov,'Style','Pushbutton','Units','normalized','Position',[0.025+woSpalten(vollSpalten) .97-((1+i-((vollSpalten-1)*30))*.028) (.9/spalten) .025],...
%             'String',{NetFlags.Covariates{i}},'HorizontalAlignment','center','FontWeight','demi','Callback',toCall,'ForegroundColor','k','FontSize',FS(8));
%           NetFlags.activeCov(i) = 0;
%         end
%       end
%     end
  else
    uicontrol(NetFlags.Stats.cov,'Style','Pushbutton','Units','normalized','Position',[0.025+woSpalten(vollSpalten) .97-((1+i-((vollSpalten-1)*30))*.028) (.9/spalten) .025],...
      'String',{NetFlags.Covariates{i}},'HorizontalAlignment','center','FontWeight','demi','Callback',toCall,'ForegroundColor','k','FontSize',FS(8));
    NetFlags.activeCov(i) = 0;
  end
end
