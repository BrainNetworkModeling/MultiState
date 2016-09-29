function STAT =se_regr(x,y)

Q = find(isnan(x) | isnan(y));
x(Q) = [];
y(Q) = [];

[x idx]=sort(x);
ytmp=y(idx);
xtmp=[x ones(length(x),1)]; %input matrix for regress function

%regression coefficients
[p,pINT,R,Rint] = regress(ytmp,xtmp);

xtmp(:,2)=[]; %delete column 2

n=length(xtmp);
xm=mean(xtmp); xsd=std(xtmp);

%standard error of regression coefficients
%Student's critical value
cv=tinv(0.975,n-2);

ystar=ytmp-R;

%Residual Variability
vres=sum(R.^2);

%regression standard error (RSE)
RSE=realsqrt(vres/(n-2));

%Confidence interval at 95% of regression
sy=RSE*realsqrt(1/n+(((xtmp-xm).^2)/((n-1)*xsd^2)));
cir=[ystar+cv*sy ystar-cv*sy];

%Confidence interval at 95% of a new observation (this is the confidence
%interval that should be used when you evaluate a new y with a new observed
%x)
sy2=realsqrt(sy.^2+RSE^2);
cir2=[ystar+cv*sy2 ystar-cv*sy2];

STAT.xtmp = xtmp;
STAT.cir = cir;
STAT.cir2 = cir2;
STAT.ystar = ystar;
