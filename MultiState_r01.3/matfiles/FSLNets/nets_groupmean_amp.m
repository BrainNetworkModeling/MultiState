%
% nets_groupmean - estimate group mean/one-group-t-test and consistency of netmats across runs/subjects
% Steve Smith, 2012-2014
%
% [Znet]      = nets_groupmean(netmats,make_figure);
% [Znet,Mnet] = nets_groupmean(netmats,make_figure);
%
% [Znet]      = nets_groupmean(netmats,make_figure,Nsubgroup);
% [Znet,Mnet] = nets_groupmean(netmats,make_figure,Nsubgroup);
%
% make_figure (0 or 1) controls whether to display the one-group-t-test group-level netmat and consistency figure
%
% Nsubgroup says that sets of Nsubgroup runs are from the same subject and should be averaged before cross-subject consistency estimation
%
% Znet is Z-stat from one-group t-test across subjects
% Mnet is mean netmat across subjects
%

function [Znet,Mnet]=nets_groupmean_amp(netmats,gofigure,varargin);

global NetFlags
NetFlags.Label = strrep(NetFlags.Label,'_',' ');

Nsubgroup=1;
if nargin==3
  Nsubgroup=varargin{1};
end

Nf=sqrt(size(netmats,2));  N=round(Nf);  Nsub=size(netmats,1);

% one-group t-test
grot=netmats; DoF=Nsub-1;
if Nsubgroup>1
  clear grot;
  for i=1:Nsub/Nsubgroup
    grot(i,:)=mean(netmats((i-1)*Nsubgroup+1:i*Nsubgroup,:));
  end
  DoF=i-1;
end

%[grotH,grotP,grotCI,grotSTATS]=ttest(grot,0);  Tnet=grotSTATS.tstat;  Tnet(isfinite(Tnet)==0)=0;
Tnet = sqrt(size(grot,1)) * mean(grot) ./ std(grot); 

Mnet=mean(grot);
Snet=std(grot);

% Znet=sign(Tnet).*(2^0.5).*erfinv(1-2.*(betainc(DoF./(DoF+abs(Tnet).^2),DoF/2,1/2)/2));
% Znet(isinf(Znet)==1)=20*sign(Znet(isinf(Znet)==1));  % very large t values would otherwise be called infinite
Znet = zeros(size(Tnet));
Znet(Tnet>0) = -norminv(tcdf(-Tnet(Tnet>0),DoF));
Znet(Tnet<0) = norminv(tcdf(Tnet(Tnet<0),DoF));

Znetd=Znet;
Mnetd=Mnet;
if N==Nf      % is netmat square....
  Znet=reshape(Znet,N,N);
  Mnet=reshape(Mnet,N,N);
end

if gofigure>0
  FIG=figure('position',[100 100 1100 400]);
  subplot(1,2,1);
  barh(Mnetd); %hold on;
%  errorbarh(Mnetd,Snet,'d'); hold off;
  if N==Nf      % is netmat square....
    Mnetd=reshape(Mnetd,N,N);
    if sum(sum(abs(Mnetd)-abs(Mnetd')))<0.00000001    % .....and symmetric
      imagesc(Mnetd,[-10 10]);
    end
  end
  set(gca,'yTick',[1:numel(NetFlags.Label)],'yTickLabel',NetFlags.Label);
    for i=1:numel(NetFlags.Label)
%      text(Mnetd(i)+0.01,i,num2str(Snet(i)));
      text(Mnetd(i)+0.01,i,['(' num2str(round(Snet(i)*1000)/1000) ')']);
    end 
  title('standard deviation of time-series per node (and it''s SD)');
  

  % scatter plot of each session's netmat vs the mean netmat
  subplot(1,2,2); 
  grot=repmat(mean(netmats),Nsub,1);
  scatter(netmats(:),grot(:));
  for i=1:numel(NetFlags.Label)
      text(round(max(netmats(:))),Mnetd(i),NetFlags.Label(i));
  end 
  title('scatter of each subjects''s std vs mean std');
end

print(FIG,'-dpng',sprintf(fullfile(pwd,'Projects','RS','Networks',spm_str_manip(NetFlags.Group,'rt'),'FIG.png')));
%print(FIG,'-append','-dpsc2',sprintf(fullfile(pwd,'Projects','RS','Networks',spm_str_manip(NetFlags.Group,'rt'),[NetFlags.project  '.ps'])));
%close(FIG);

end