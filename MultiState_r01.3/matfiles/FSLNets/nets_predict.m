%
% nets_predict - predict data using multiple regression, with feature (design matrix column) pre-selction, using LOO and permutation testing
% Steve Smith, Diego Vidaurre, Ludo Griffanti, Tom Nichols, Anderson Winkler
% FMRIB Oxford, 2013-2014
%
% [pval] = nets_predict(Y,X,Nfeatures,Nperm,show_scatter,);
% [pval] = nets_predict(Y,X,Nfeatures,Nperm,show_scatter,correlation_structure);
% [pval] = nets_predict(Y,X,Nfeatures,Nperm,show_scatter,correlation_structure,Permutations);
% [pval,predictedY] = nets_predict()
%
% INPUTS
% Y - data vector (samples X 1)
% X - design matrix (samples X features)
% Nfeatures - number of strongest features to keep (set to 0 to use all features)
% Nperm - number of permutations (set to 0 to skip permutation testing)
% show_scatter - set to 1 to show a scatter plot of predicted_Y vs Y
% correlation_structure (optional) - an (Nsamples X Nsamples) matrix with integer dependency labels (e.g., family structure)
% Permutations (optional but must also have correlation_structure) - pre-created set of permutations
%
% OUTPUTS
% pval - permutation-based p-value, if permutation is run; otherwise, correlation-based p-value
% predictedY (optional) - predicted data vector
%

function [pval,predictedY,rval] = nets_predict(Yin,Xin,Nfeatures,Nperm,show_scatter,varargin);

if (Nperm<2),  Nperm=1;  end;
if (Nfeatures<1),  Nfeatures=0;  end;
if (Nfeatures>=size(Xin,2)),  Nfeatures=0;  end;
N=size(Yin,1);  YinORIG=Yin;
cs=[]; if (nargin>5)
  cs=varargin{1};
  [allcs(:,2),allcs(:,1)]=ind2sub([length(cs) length(cs)],find(cs>0));   % find all correlated samples for every sample
  [grotMZi(:,2),grotMZi(:,1)]=ind2sub([length(cs) length(cs)],find(tril(cs,1)==1));
  [grotDZi(:,2),grotDZi(:,1)]=ind2sub([length(cs) length(cs)],find(tril(cs,1)==2));
  %[grotSIi(:,2),grotSIi(:,1)]=ind2sub([length(cs) length(cs)],find(tril(cs,1)==3));
end
PrePerms=0; if (nargin>6)
  Permutations=varargin{2};
  PrePerms=1;
  Nperm=size(Permutations,2);
end

for perm=1:Nperm
  if (perm>1)
    perm
    if (length(cs)==0)            % simple full permutation with no correlation structure
      Yin=Yin(randperm(N));
    elseif (PrePerms==0)          % complex permutation, taking into account correlation structure
      PERM=zeros(1,N);
      perm1=randperm(size(grotMZi,1));
      for ipe=1:length(perm1)
        if rand<0.5, wt=[1 2]; else, wt=[2 1]; end;
        PERM(grotMZi(ipe,1))=grotMZi(perm1(ipe),wt(1));
        PERM(grotMZi(ipe,2))=grotMZi(perm1(ipe),wt(2));
      end
      perm1=randperm(size(grotDZi,1));
      for ipe=1:length(perm1)
        if rand<0.5, wt=[1 2]; else, wt=[2 1]; end;
        PERM(grotDZi(ipe,1))=grotDZi(perm1(ipe),wt(1));
        PERM(grotDZi(ipe,2))=grotDZi(perm1(ipe),wt(2));
      end
      from=find(PERM==0);  pto=randperm(length(from));  to=from(pto);  PERM(from)=to;
      Yin=YinORIG(PERM);
    else                   % pre-supplied permutation
      Yin=YinORIG(Permutations(:,perm));  % or maybe it should be the other way round.....?
    end
  end

  grotDONE=zeros(1,N);
  for j=1:N  % LOO loop
   if (grotDONE(j)==0)
    J=j;
    if (length(cs)>0)  % leave out all samples related (according to cs) to the one in question
      if size(find(allcs(:,1)==j),1)>0
        J=[j allcs(find(allcs(:,1)==j),2)'];
      end
    end
    grotDONE(J)=1;
    ji=setdiff(1:N,J);  X=nets_normalise(Xin(ji,:));  my=mean(Yin(ji)); Y=Yin(ji)-my;
    if (Nfeatures>0)
      [~,groti]=sort(abs(Y'*X));  groti=groti(end-Nfeatures+1:end);  X=X(:,groti);
    else
      groti=1:size(Xin,2);
    end
    mx=mean(Xin(ji,groti));  sx=std(Xin(ji,groti));  grotbeta=pinv(X)*Y;
    for jj=J
      predictedYp(jj)=((Xin(jj,groti) - mx) ./ sx)  * grotbeta + my;
    end
   end
  end

  grotc=(Yin-predictedYp'); grotperms(perm)=sqrt(mean(grotc.*grotc));
  if perm==1
    predictedY=predictedYp';
    if show_scatter
      figure;  scatter(Yin,predictedY);
    end
  end
end

[grotr,pval] = corrcoef(YinORIG,predictedY);  pval=pval(1,2);  rval=grotr(1,2)

if (Nperm>1)
  pval = sum(grotperms<=grotperms(1)) / Nperm;
end

