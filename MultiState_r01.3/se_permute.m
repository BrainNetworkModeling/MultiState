function p = se_permute(xdata,groups,type)

repeats = 25001;


try
    test = ones(1,10); parfor i=1:10; test(i) = 0; end
    useparellel = 1;
catch
    useparellel = 0;
end

        
if type == 1
  benchmark = nanmedian(xdata(groups>.5))-nanmedian(xdata(groups<.5));
  A = 0;
  B = 0;
  q = numel(xdata);
  if useparellel
      parfor i=1:repeats
          data = xdata(randperm(q));
          A = A + double((nanmedian(data(groups>.5))-nanmedian(data(groups<.5)))>benchmark);
      end
  else
      for i=1:repeats
          data = xdata(randperm(q));
          A = A + double((nanmedian(data(groups>.5))-nanmedian(data(groups<.5)))>benchmark);
      end
  end
  
else
  benchmark = nanmean(xdata(groups>.5))-nanmean(xdata(groups<.5));
  A = 0;
  B = 0;
  q = numel(xdata);
  if useparellel
      parfor i=1:repeats
          data = xdata(randperm(q));
          A = A + double((nanmean(data(groups>.5))-nanmean(data(groups<.5)))>benchmark);
      end
  else
      for i=1:repeats
          data = xdata(randperm(q));
          A = A + double((nanmean(data(groups>.5))-nanmean(data(groups<.5)))>benchmark);
      end
  end
end

p = min(A/repeats,(repeats-A)/repeats);