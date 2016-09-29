Vi = spm_vol(spm_select(Inf,'image','Select image files'));


titlestr = spm_input('Project Name','+1','s');

for i=1:numel(Vi)
  dat = spm_read_vols(Vi(i));
  ind = find(dat>0);
  [X Y Z] = ind2sub(Vi(i).dim,ind);
  nXYZ{i} = [X Y Z]';
  clear X Y dat ind;
end


fid = fopen(fullfile(pwd,'VOIs',[titlestr '_ImageList.txt']),'wt');
for i=1:numel(Vi)
  fg = spm_figure('GetWin','Graphics'); spm_figure('Clear','Graphics'); spm_orthviews('Reset');
  try
    spm_orthviews('Image', spm_vol('Anatomy_V20alpha_MNI.nii'), [0.0 0.22 1 .8]);
  catch
    spm_orthviews('Image', spm_vol('ColinNIFTI_MNI.nii'), [0.0 0.22 1 .8]);
  end
  
  for ii=1:numel(Vi)
    if ii~=i
      spm_orthviews('addcolouredblobs',1,nXYZ{ii},1,Vi(ii).mat,[1 0 0]);
    end
  end
  mXYZ   = mean(Vi(i).mat * [nXYZ{i}; ones(1,size(nXYZ{i},2))],2); mXYZ   = mXYZ(1:3)';
  spm_orthviews('addcolouredblobs',1,nXYZ{i},1,Vi(i).mat,[0 1 0]);
  spm_orthviews('reposition',mXYZ)
  ROIname = spm_input('VOI name [0 skips]',1,'s',[spm_str_manip(Vi(i).fname,'rt') ]);
  fprintf(fid,'%s\t%s\n',Vi(i).fname,ROIname);
end
fclose(fid);


OverGUI