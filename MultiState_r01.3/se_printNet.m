function se_printNet

global NetFlags

if strcmp(NetFlags.Stats.Results{1}(1:6),'No sig')
  spm('alert*','No results to print')
else
  fname = spm_input('File name','+1','s',[spm_str_manip(NetFlags.Group,'rt') '_' NetFlags.project '_' NetFlags.Stats.Type]);
  for i=1:numel(NetFlags.Stats.Results)
    se_ShowEffect(i)
    F = findall(allchild(0),'Flat','Tag','Show Effect');
    set(0,'CurrentFigure',F);
    figure(F);
    print(F,[fname '.ps'],'-append','-dpsc2')
  end
end
