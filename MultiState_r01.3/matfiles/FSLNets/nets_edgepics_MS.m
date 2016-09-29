%
% nets_edgepics(ts,group_maps,meannetmat,net,showN);
% Steve Smith - 2013-2014
%
% show the strongest netmat elements from "net" (the number shown is controlled by showN)
% 
% meannetmat is typically the group-mean or group-mean-Z netmat, used just to generate the coloured bars connecting the two nodes in each edge
%
% ts is the node timeseries structure; group_maps is a string pointing to the thumbnail images folder
%

function [netmat] = nets_edgepics(ts,group_maps,meannetmat,net,showN,des)    %%%% show a snapshot of the kept and rejected components

global NetFlags;

Label = strrep(NetFlags.Label,'_',' ');
netORIG=net;
net=abs(net);
net(tril(ones(size(net,1)))==1)=0;   % get info on upper tri only
meannetmat=max(min(meannetmat/prctile(abs(meannetmat(:)),90),1),-1);
grott=sprintf('%s.png',tempname);

[grotA,grotB,grotC]=fileparts(group_maps); if size(grotA,1)==0, grotA='.'; end; group_maps=sprintf('%s/%s.sum',grotA,grotB);

[yy,ii]=sort(net(:),'descend');   % find strongest showN edges

xf=ceil(1.15*sqrt(showN)); yf=ceil(showN/xf);   % dimensions of display tiling

XX=1000; YY= 500; % total image size

bd=15; bdx=bd/XX; bdy=bd/YY; % borders size
el=10/XX; % edge lengths
th=20/YY; % title height
isx=(1-(xf+1)*bdx-xf*el)/(2*xf);  isy=(1-(yf+1)*bdy-yf*th)/yf;

YY=(109/91) * (isx*XX) / isy;
FIG=figure('Position',[10 10 XX YY]); 

for iii=1:showN
  xxx=ceil(ii(iii)/size(net,1));  yyy=ii(iii)-((xxx-1)*size(net,1));
  yfi=ceil(iii/xf); xfi=iii-((yfi-1)*xf); yfi=1+yf-yfi;   [xxx yyy xfi yfi];

  system(sprintf('slices_summary %s %s %s',group_maps,grott,num2str(ts([xxx])-1)));  picgood=imread(grott);
  subplot('position',[xfi*bdx+(xfi-1)*(2*isx+el) yfi*bdy+(yfi-1)*(isy+th) isx isy]);
  imagesc(picgood); axis off;

  system(sprintf('slices_summary %s %s %s',group_maps,grott,num2str(ts([yyy])-1)));  picgood=imread(grott);
  subplot('position',[xfi*bdx+(xfi-1)*(2*isx+el)+isx+el yfi*bdy+(yfi-1)*(isy+th) isx isy]);
  imagesc(picgood); axis off; %axis equal;

  grot=meannetmat(yyy,xxx); eh=abs(grot)*0.5*isy; [grot eh];
  subplot('position',[xfi*bdx+(xfi-1)*(2*isx+el)+isx yfi*bdy+(yfi-1)*(isy+th)+isy/2-eh/2 el eh]);
  imagesc(sign(grot)-0.5,[-.8 .8]); axis off;

  subplot('position',[xfi*bdx+(xfi-1)*(2*isx+el) yfi*bdy+(yfi-1)*(isy+th)+isy 2*isx+el th]); axis off;
  title(sprintf('%s <-> %s   Z=%.2f',Label{xxx},Label{yyy},netORIG(yyy,xxx)),'Position',[0.5 0.5]);
end

set(gcf,'PaperPositionMode','auto','PaperOrientation','landscape');


print(FIG,'-append','-dpsc2',sprintf(fullfile(pwd,'Projects','RS','Networks',spm_str_manip(NetFlags.Group,'rt'),[NetFlags.project '_' des],['clustering.ps'])));
%close(FIG);

end