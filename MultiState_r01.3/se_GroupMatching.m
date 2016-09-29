function [exPAT, exCON] = se_GroupMatching(isPat,nPAT,nCON,Gender,Features)


critical = 0.2;

if numel(unique(Gender([nPAT nCON])))==1
    doGender = false;
elseif numel(unique(Gender([nPAT nCON])))==2
    doGender = true;
else
    error('More than 2 categories of gender - are you sure?');
end

todo = size(Features,2);
ps = zeros(1,todo+1);
if doGender
    [~,~,ps(1)] = crosstab(isPat([nPAT,nCON]), Gender([nPAT nCON]));
else
    ps(1) = 1;
end

for i=1:todo; [~, ps(i+1)] = ttest2(Features(nPAT,i),Features(nCON,i));   end
if all(ps>critical); exPAT = [];  exCON = []; fprintf(1,'%s\n','OK'); return; end

[B I] = sort(ps(2:end),'ascend');
Features = Features(:,I);

exPAT  = []; np = numel(nPAT);
exCON  = []; nc = numel(nCON);

for i=1:max(floor(min(np,nc)/4),10)
    np = numel(nPAT); nc = numel(nCON);
    if min(np,nc)<12
        break
    end
    
    ps = zeros(1,todo+1); if doGender ; [~,~,ps(1)] = crosstab(isPat([nPAT,nCON]), Gender([nPAT nCON])); else, ps(1) = 1; end
    for ii=1:todo; [~, ps(ii+1)] = ttest2(Features(nPAT,ii),Features(nCON,ii));   end
    if min([ps])>.005
        break
    end
    
    outoP = zeros(np,todo);
    outoC = zeros(nc,todo);
    for ii=1:todo
        if ps(ii+1)<0.005
            A = pdist2(Features(nPAT,ii),Features(nCON,ii));
            if max(mean(A,2))>max(mean(A,1))
                [B I] = sort(mean(A,2)); outoP(I,ii) = 1:numel(nPAT);
            else
                [B I] = sort(mean(A,1)); outoC(I,ii) = 1:numel(nCON);
            end
        end
    end
    if max(sum(outoP,2))>max(sum(outoC,2))
        [B I] = sort(sum(outoP,2),'descend'); exPAT = [exPAT nPAT(I(1))]; nPAT(I(1)) = [];
    else
        [B I] = sort(sum(outoC,2),'descend'); exCON = [exCON nCON(I(1))]; nCON(I(1)) = [];
    end
end


ps = zeros(1,todo+1); if doGender ; [~,~,ps(1)] = crosstab(isPat([nPAT,nCON]), Gender([nPAT nCON])); else, ps(1) = 1; end
for i=1:todo; [~, ps(i+1)] = ttest2(Features(nPAT,i),Features(nCON,i));   end

if all(ps>critical); exPAT = [];  exCON = []; fprintf(1,'%s\n','OK'); return; end

[B I] = sort(ps(2:end),'ascend');
Features = Features(:,I);



repeats = 1E6; tic
np = numel(nPAT);
nc = numel(nCON);


keepPAT = {};
keepCON = {};
keepAnz = [];
keepPs  = [];
keepPs2  = [];

tic
for offset = 0:7;
    
    for depth = 3:6;
        for iteration = 1:10
            minanz = max(min(numel(nPAT),numel(nCON))-depth-offset,10);
            pps = floor(rand(repeats*3,1)*(numel(nPAT)-minanz+1))+minanz-offset;
            ccs = floor(rand(repeats*3,1)*(numel(nCON)-minanz+1))+minanz-offset;
            Q = (ccs./pps)>1.2 | (pps./ccs)>1.2  | abs(pps-ccs)>5;
            pps(Q) = []; pps = pps(1:min(repeats,sum(Q==0)));
            ccs(Q) = []; ccs = ccs(1:min(repeats,sum(Q==0)));
            
            repeats = min(repeats,min(numel(pps),numel(ccs)));
            pps = pps(1:repeats);
            ccs = ccs(1:repeats);
            
            
            % Randomly select subjects but remove duplicates
            
            mp = max(pps); allNowP = zeros(repeats,mp);
            mc = max(ccs); allNowC = zeros(repeats,mc);
            
            parfor rep = 1:repeats
                allNowP(rep,:) = [randsample(nPAT,pps(rep)) zeros(1,mp-pps(rep))];
                allNowC(rep,:) = [randsample(nCON,ccs(rep)) zeros(1,mc-ccs(rep))];
            end
            clear  ccs pps
            
            allNowP = sort(allNowP,2,'descend'); allNowC = sort(allNowC,2,'descend');
            allNow = unique([allNowP allNowC],'rows');
            
            allNowP = uint16(allNow(:,1:mp)); allNowC = uint16(allNow(:,mp+1:end));
            clear allNow
            repeats = size(allNowP,1);
            
            
            % Evaluate random selections
            emptyp  = zeros(1,todo+1);
            usePAT  = cell(repeats,1);
            useCON  = cell(repeats,1);
            allPs   = zeros(repeats,todo+1);
            total   = zeros(repeats,1);
            
            parfor rep = 1:repeats
                nowP = allNowP(rep,:); nowP(nowP==0) = [];
                nowC = allNowC(rep,:); nowC(nowC==0) = [];
                nP = emptyp;
                for i=1:todo
                    [~, nP(i+1)] = ttest2(Features(nowP,i),Features(nowC,i));
                    if nP(i+1)<critical,
                        break;
                    end
                end
                if all(nP(2:end)>critical)
                    if doGender
                        [~,~,nP(1)] = crosstab(isPat([nowP,nowC]), Gender([nowP,nowC]));
                    else
                        nP(1) = 1;
                    end
                end
                if all(nP>critical) & min(numel(nowP),numel(nowC))>9
                    total(rep) = numel([nowP nowC])+min(numel(nowP),numel(nowC));
                    usePAT{rep} = nowP;
                    useCON{rep} = nowC;
                    allPs(rep,:) = nP;
                end
            end
            
            if max(total)>0
                Q = (total>=max(total)-1);
                Q = find(Q);
                for xi=1:numel(Q)
                    keepPAT{end+1}   = usePAT{Q(xi)};
                    keepCON{end+1}   = useCON{Q(xi)};
                    keepAnz(end+1,:) = [numel(usePAT{Q(xi)}) numel(useCON{Q(xi)})];
                    keepPs(end+1)    = min(allPs(Q(xi),:));
                    keepPs2(end+1)   = mean(allPs(Q(xi),:));
                end
            end
        end
        fprintf(1,'%s\n',['10 iterations @ 1E6 of depth ' int2str(depth) ' / offset ' int2str(offset)  ' finished after ' int2str(toc/60) ' minutes']);
    end
    if numel(keepPs)>10 & offset>1
        break
    end
end


fprintf(1,'%s\t',['       ' int2str(toc) 's']);
if numel(keepPs)==0
    exPAT = [exPAT nPAT];
    exCON = [exCON nCON];
    fprintf(1,'%s\n',[' ALL excluded']);
else
    total = min(keepAnz,[],2) + max(keepAnz,[],2)*.33;
    Q = find(total>=max(total)-.5);
    Q = Q(find(keepPs(Q)==max(keepPs(Q))));
    if numel(Q)>0; Q = Q(total(Q)==max(total(Q))); end
    if numel(Q)>0; Q = Q(find(keepPs2(Q)==max(keepPs2(Q)))); end
    Q = Q(1);
    
    try
        for i=1:np; if ~any(keepPAT{Q}==nPAT(i)); exPAT = [exPAT nPAT(i)]; end; end
        for i=1:nc; if ~any(keepCON{Q}==nCON(i)); exCON = [exCON nCON(i)]; end; end
        fprintf(1,'%s\n',[int2str(numel(exPAT)) ' / ' int2str(numel(exCON)) ' excluded']);
    catch
        np
    end
end

end
