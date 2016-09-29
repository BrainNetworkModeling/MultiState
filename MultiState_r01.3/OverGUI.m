clc, clear,

clear global NetFlags GLMFlags

fg = spm_figure('getWin','Graphics');
spm_figure('Clear','Graphics');

set(gcf,'DefaultUicontrolFontSize',spm('FontSizes',get(gcf,'DefaultUicontrolFontSize')));
WS = spm('WinScale');


uicontrol(fg,'Style','Text','Units','normalized','Position',[.05 .9 .9 .08],'String','The multi-state box',...
  'FontSize',spm('FontSizes',24),'FontWeight','bold','BackgroundColor',[1 1 1]);

uicontrol(fg,'Style','PushButton','Units','normalized','Position',[.05 .85 .9 .05],'Callback','se_LookupFromExcel;',...
  'String','Create Lookup-file from Excel spreadsheet','FontSize',spm('FontSizes',14));

uicontrol(fg,'Style','PushButton','Units','normalized','Position',[.05 .75 .9 .05],'Callback','se_SetupVOIsAsText;',...
  'String','Setup VOI coordinate textfile','FontSize',spm('FontSizes',14));

uicontrol(fg,'Style','PushButton','Units','normalized','Position',[.05 .7 .9 .05],'Callback','se_Visualize;',...
  'String','Visualize VOI coordinate locations','FontSize',spm('FontSizes',14));


uicontrol(fg,'Style','PushButton','Units','normalized','Position',[.05 .625 .9 .05],'Callback','se_splitImage;',...
  'String','Create individual binary VOIs from image','FontSize',spm('FontSizes',14));

uicontrol(fg,'Style','PushButton','Units','normalized','Position',[.05 .575 .9 .05],'Callback','se_SetupImageVOIs;',...
  'String','Setup project from VOI images','FontSize',spm('FontSizes',14));

uicontrol(fg,'Style','PushButton','Units','normalized','Position',[.05 .5 .9 .05],'Callback','se_Network;',...
  'String','Network RS connectivity analysis','FontSize',spm('FontSizes',14));

uicontrol(fg,'Style','PushButton','Units','normalized','Position',[.05 .45 .9 .05],'Callback','se_GLM;',...
  'String','Whole-brain RS connectivity analysis','FontSize',spm('FontSizes',14));


uicontrol(fg,'Style','PushButton','Units','normalized','Position',[.05 .375 .9 .05],'Callback','se_Volumes;',...
  'String','Regional VBM analysis','FontSize',spm('FontSizes',14));

uicontrol(fg,'Style','PushButton','Units','normalized','Position',[.05 .325 .9 .05],'Callback','se_GLMVBM2;',...
  'String','Whole-brain VBM incl Structural Covariance','FontSize',spm('FontSizes',14));


uicontrol(fg,'Style','PushButton','Units','normalized','Position',[.05 .250 .9 .05],'Callback','se_TBSS;',...
  'String','TBSS whole-brain FA Analysis','FontSize',spm('FontSizes',14));

