clc

if ~(exist('VOI_Locations','dir')==7)
  mkdir('VOI_Locations');
end

TxtName       = spm_select(1,'any',['Select coordinate list'],[],fullfile(pwd,'VOIs'),'.txt',1);
[X Y Z Label]     = textread(TxtName,'%f %f %f %s');
project = strrep(spm_str_manip(TxtName,'rt'),'_VOIs','');

load ./MaskenEtc/MNI152GM

plotXYZ = [];
T       = [];

sXYZ = VM(1).mat * [maskXYZ; ones(1,size(maskXYZ,2))];

for ana = 1:numel(X)
  
  D    = sqrt((sXYZ(1,:)-X(ana)).^2+(sXYZ(2,:)-Y(ana)).^2+(sXYZ(3,:)-Z(ana)).^2);
  plotXYZ = [plotXYZ maskXYZ(:,D<5)];
  T       = [T ones(1,sum(D<5))*ana];
end

dat = accumarray(plotXYZ',1,VM.dim);

Vo = VM;
Vo = rmfield(Vo,'pinfo');

Vo.fname = fullfile(pwd,'VOI_Locations',[project '.nii']);
spm_write_vol(Vo,dat);

se_render_imageCol(fullfile(pwd,'VOI_Locations',[project '.nii']),0,0,[1 0 0])
print('-dpng',fullfile(pwd,'VOI_Locations',[project '_3D.png']))
ImageCut(fullfile(pwd,'VOI_Locations',[project '_3D.png']),'X')
delete(fullfile(pwd,'VOI_Locations',[project '_3D.png']));



plotXYZ = [plotXYZ; ones(1,size(plotXYZ,2))];



setp = [...
  1 2 2 2 3 3 3 3 3 4 4 4 4 4 4 4 5 5 5 5 6 6 6 6  
  1 1 2 2 2 2 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 4 4];
  

setp = setp(:,max(T));

fg = spm_figure('GetWin','Graphics'); spm_figure('Clear','Graphics'); spm_orthviews('Reset');
spm_orthviews('Image', spm_vol('Anatomy_V20alpha_MNI.nii'), [0.0 0.22 1 .8]);
for i=1:max(T)
    spm_orthviews('addcolouredblobs',1,plotXYZ(1:3,T==i),1,Vo.mat,[1 0 0]);
    mXYZ = median([Vo.mat * double(plotXYZ(:,T==i))]');
    spm_orthviews('reposition',mXYZ(1:3));
    print('-dpng',fullfile(pwd,'VOI_Locations',['TMP_' int2str(i) '.png']));
    se_orthviews('rmblobs',1);
end

for i=1:max(T)
    A = imread(fullfile(pwd,'VOI_Locations',['TMP_' int2str(i) '.png']),'png');
    delete(fullfile(pwd,'VOI_Locations',['TMP_' int2str(i) '.png']));
    figure(90); subplot(setp(1),setp(2),i); imshow(A(100:560,1:460,:)); title(strrep(Label{i},'_',' '));
end
was = get(gcf,'Position'); set(gcf,'Position',get(0,'ScreenSize')); 
print('-dpng',fullfile(pwd,'VOI_Locations',[project '_Slices.png']));
delete(gcf);
spm_figure('Clear','Graphics'); 

OverGUI
