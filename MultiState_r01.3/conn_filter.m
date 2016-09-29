function y=conn_filter(rt,filter,x,option);

if nargin<4, option='full'; end
y=fft(x,[],1);
f=(0:size(x,1)-1);
f=min(f,size(x,1)-f);
switch(lower(option))
	case 'full',
		idx=find(f<filter(1)*(rt*size(x,1))|f>filter(2)*(rt*size(x,1)));
		idx=idx(idx>1);
		y(idx,:)=0;
		y=real(ifft(y,[],1));
	case 'partial',
		idx=find(f>=filter(1)*(rt*size(x,1))&f<=filter(2)*(rt*size(x,1)));
		if ~any(idx==1), idx=[1,idx]; end
		y=y(idx,:);
end

