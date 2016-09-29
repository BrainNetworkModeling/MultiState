clc, clear

gname = spm_select(Inf,'any','Select lookup speadsheets','',fullfile(pwd,'Phenotypical'),'.xls');
[A B] = xlsread(gname,1,'','basic');

who = spm_input('Include','!+1','b','All|RS/VBM|3Modal',[1 2 3],1,'Select mode');
wat = spm_input('With Cov ->','!+1','b','Complete|All',[1 2],1,'Select mode');
match = spm_input('Match groups','!+1','b','Movements/Age/Sex|Age/Sex|None',[1 2 0],1,'Select mode');
% out = spm_input('Covariate outlier exclusion','!+1','b','> 2 std | > 3 std | None',[2 3 0],1,'Select mode');
% if out > 0
%     outrow = spm_input('1st XLS-row for outliers','!+1');
% end
fid = fopen('_path2data','r'); datadir = fgetl(fid); fclose(fid);

if A(1,1)~=99 % & isempty(B{2,1})
    Q = find(~isnan(A(:,1)));
    for i=1:numel(Q)
        B{Q(i),2} = int2str(A(Q(i),1));
    end
    A = A(:,2:end);
end

Covariates = {B{1,[4:end]}};

for i=2:size(A,1)
    xsub{i-1} = B{i,2};
    xSubDir{i-1} = B{i,1};
    if exist(fullfile(datadir,'DATA',xSubDir{i-1},xsub{i-1}))==0
        if exist(fullfile(datadir,'DATA',xSubDir{i-1},['0' xsub{i-1}]))>0
            xsub{i-1} = ['0' xsub{i-1}];
        elseif  exist(fullfile(datadir,'DATA',xSubDir{i-1},['00' xsub{i-1}]))>0
            xsub{i-1} = ['00' xsub{i-1}];
        end
    end
end

hasRS  = [];
hasRSfix= [];
hasVBM = [];
hasDTI = [];

fprintf(1,'%s\n',['Checking imaging data of ' int2str(numel(xsub)) ' subjects ...']);
cnt = 1;
for i=1:numel(xsub)
    
    if exist(fullfile(datadir,'DATA',xSubDir{i},xsub{i},'RS',['MNI152GM_' xsub{i} '.nii']));
        hasRS(i) = 1;
        tmp = load(fullfile(datadir,'DATA',xSubDir{i},xsub{i},'RS',['Counfounds_' xsub{i} '.mat']));
        scans(i) = numel(tmp.Vi);
    else
        hasRS(i) = 0;
        scans(i) = NaN;
    end
    
    if exist(fullfile(datadir,'DATA',xSubDir{i},xsub{i},'RS_fix',['MNI152GM_' xsub{i} '.nii']));
        hasRSfix(i) = 1;
    else
        hasRSfix(i) = 0;
    end
    
    if exist(fullfile(datadir,'DATA',xSubDir{i},xsub{i},'3D',['m0wrp1' xsub{i} '.nii']));
        try
            hasVBM(i) = 1;
            tbv(i,[1 2 3]) = load(fullfile(datadir,'DATA',xSubDir{i},xsub{i},'3D',['p' xsub{i} '_seg8.txt']));
        catch
            hasVBM(i) = 0;
            tbv(i,[1 2 3]) = nan(1,3);
        end
    else
        hasVBM(i) = 0;
        tbv(i,[1 2 3]) = nan(1,3);
    end
    
    if exist(fullfile(datadir,'DATA',xSubDir{i},xsub{i},'DWI','reg3G','Mean3G_warp2FA.nii.gz')) && ...
            exist(fullfile(datadir,'DATA',xSubDir{i},xsub{i},'DWI',[xsub{i} '.bedpostX'],'mean_f1samples.nii.gz'));
        hasDTI(i) = 1;
    else
        hasDTI(i) = 0;
    end
    
    sub{cnt} = xsub{i};
    SubDir{cnt} = xSubDir{i};
    
    Cov(cnt,:) = A(i+1,2:end);
    isPat(cnt) = A(i+1,1);
    cnt = cnt+1;
    if rem(i,10)==0 
        fprintf(1,'%s','.');
    end
end

if sum(hasRS)+sum(hasVBM)+sum(hasDTI)==0
    error('No imaging data found, check path2date & your phenotypical');
end

%%%     identify scanning sites
if any(strcmp(Covariates,'Site'))
    si = find(strcmpi(Covariates,'Site')==1);
    sites = unique(Cov(:,si));
    scanner = Cov(:,si);
else
    scanner(1:numel(sub),1) = 1;
    sites = 1;
end
   
%%%     exclude subjects with less than 95% time points of site majority
for i = 1:numel(sites)
    minTP(i) = prctile(scans(scanner==sites(i)),10)*.95;
    excl{i} = find(scans(scanner==sites(i))<minTP(i));
    hasRS(excl{i})=0;
    if numel(excl{i})>0
        spm('alert"',{['Site ' int2str(i) ': ' int2str(numel(excl{i}))  ' subjects with less than ' ...
            int2str(minTP(i)) ' time points excluded']});
    end
end

if who==1
    use = sum([hasRS; hasVBM; hasDTI])>0;
    suffix = 'all';
elseif who == 2
    use = (hasRS+hasVBM)==2;
    suffix = 'RSandVBM';
else
    use = sum([hasRS; hasVBM; hasDTI])==3;
    suffix = 'RSandVBMandDWI';
end

adfix = '';
if any(isPat==1)
    if match == 1
        adfix = '_mov_matched';
    elseif match == 2
        adfix = '_matched';
    end
else
	  match = 0;    
end

if wat == 1
    use = use & sum(isnan(Cov(:,sum(isnan(Cov))<size(Cov,1))),2)'==0;
end


% 
% if out>0
%     outCov(1:numel(sub),1:numel(Covariates)-outrow)=0;
%     for i= outrow:numel(Covariates)
%         outCovariates{i+1-outrow} = Covariates{i+1-outrow};
%         meanC(i+1-outrow) = nanmean(Cov(isPat==0,i));
%         stdC(i+1-outrow)  = nanstd(Cov(isPat==0,i));
%         meanP(i+1-outrow) = nanmean(Cov(isPat==1,i));
%         stdP(i+1-outrow)  = nanstd(Cov(isPat==1,i));
%         outCov(find(Cov(:,i) > meanC(i+1-outrow)+out*stdC(i+1-outrow) & isPat(:)==0) ,i) = 1;
%         outCov(find(Cov(:,i) < meanC(i+1-outrow)-out*stdC(i+1-outrow) & isPat(:)==0) ,i) = 1;
%         outCov(find(Cov(:,i) > meanP(i+1-outrow)+out*stdP(i+1-outrow) & isPat(:)==1) ,i) = 1;
%         outCov(find(Cov(:,i) < meanP(i+1-outrow)-out*stdP(i+1-outrow) & isPat(:)==1) ,i) = 1;
%     end
% end


%%% define filter here if needed
% use = use & Cov(:,1)'<=75 & Cov(:,1)'>=20 & ~(Cov(:,3)'>=9) & isPat<1;

sub = sub(use);
SubDir = SubDir(use);
Cov = Cov(use,:);
hasRS  = hasRS(use);
hasRSfix  = hasRSfix(use);
hasVBM = hasVBM(use);
hasDTI = hasDTI(use);
isPat  = isPat(use);
tbv    = tbv(use,:);
scanner = scanner(use);


RMS   = nan(1,numel(sub));
FD    = nan(1,numel(sub));
DVARS = nan(1,numel(sub));
GMmean = nan(1,numel(sub));

parfor i=1:numel(sub)
    
    if hasRS(i)>0
        
        info = load(fullfile(datadir,'DATA',SubDir{i},sub{i},'RS',['Counfounds_' sub{i} '.mat']));
        rp = info.rp;
        
        drp = [zeros(1,6); diff(rp)];
        move = sqrt(sum(drp(:,1:3).^2,2));
        
        euler = zeros(size(rp,1),1);
        for ii=1:size(rp,1)
            euler(ii) = 50*acos((cos(drp(ii,4))*cos(drp(ii,5)) + cos(drp(ii,4))*cos(drp(ii,6)) + ...
                cos(drp(ii,5))*cos(drp(ii,6)) + sin(drp(ii,4))*sin(drp(ii,6))*sin(drp(ii,5)) - 1)/2);
        end
        
        RMS(i)  = sqrt(nanmean([move+euler].^2));
        tmp   = sum(abs([drp(:,1:3) drp(:,4:6)*50]),2);
        FD(i)   = sqrt(nanmean(tmp.^2));
        
        dat = spm_read_vols(spm_vol(fullfile(datadir,'DATA',SubDir{i},sub{i},'RS',['MNI152GM_' sub{i} '.nii'])));
        dat = (diff(dat)./dat(1:end-1,:))*100;
        dat(abs(dat)>10) = NaN;
        DVARS(i) = sqrt(nanmean(nanmean(dat.^2,2)));
        
    end
    
    if hasVBM(i)>0
        dat2 = spm_read_vols(spm_vol(fullfile(datadir,'DATA',SubDir{i},sub{i},'3D',['m0wrp1' sub{i} '.nii'])));
        GMmean(i) = mean(dat2(dat2>0));
    end
    
    fprintf(1,'%s\n',[int2str(i) ' / ' int2str(numel(sub))]);

end

if match >= 1
    
    ag = find(strcmpi(Covariates,'Age')==1);
    ge = find(cellfun('length',strfind(lower(Covariates),'sex'))+ ...
        cellfun('length',strfind(lower(Covariates),'gender')));
    ge = ge(1);
    
    if who == 1
        for site = 1:numel(sites)
            nPAT = find(isPat==1 & scanner' == sites(site) & hasVBM==1);
            nCON = find(isPat==0 & scanner' == sites(site) & hasVBM==1);
            if min(numel(nPAT),numel(nCON))>9
                str = repmat(' ',1,25); str2 = [SubDir{nPAT(1)} ' - VBM:']; str(1:numel(str2)) = str2;
                fprintf(1,'%s',str);
               [exPAT, exCON] = se_GroupMatching(isPat,nPAT,nCON,Cov(:,ge),Cov(:,ag));
               hasVBM([exPAT exCON]) = 0;
            else
                hasVBM([nPAT nCON]) = 0;
            end
            
            nPAT = find(isPat==1 & scanner' == sites(site) & hasRS==1);
            nCON = find(isPat==0 & scanner' == sites(site) & hasRS==1);
            if min(numel(nPAT),numel(nCON))>9
                str = repmat(' ',1,25); str2 = [SubDir{nPAT(1)} ' - RS:']; str(1:numel(str2)) = str2;
                fprintf(1,'%s',str);
                if match == 1
                   [exPAT, exCON] = se_GroupMatching(isPat,nPAT,nCON,Cov(:,ge),[Cov(:,ag) DVARS' FD' RMS']);
                    hasRS([exPAT exCON]) = 0;
                    hasRSfix([exPAT exCON]) = 0;
                else min(numel(nPAT),numel(nCON))>9;
                   [exPAT, exCON] = se_GroupMatching(isPat,nPAT,nCON,Cov(:,ge),Cov(:,ag));
                    hasRS([exPAT exCON]) = 0;
                    hasRSfix([exPAT exCON]) = 0;
                end
            else
                hasRS([exPAT exCON]) = 0;
                hasRSfix([nPAT nCON]) = 0;
            end
        end
        fprintf(1,'\n');

    else
        for site = 1:numel(sites)
            nPAT = find(isPat==1 & scanner' == sites(site) & hasRS==1);
            nCON = find(isPat==0 & scanner' == sites(site) & hasRS==1);
            if min(numel(nPAT),numel(nCON))>9
                str = repmat(' ',1,25); str2 = [SubDir{nPAT(1)} ' within site matching: - RS & VBM']; str(1:numel(str2)) = str2;
                fprintf(1,'%s',str);
                if match == 1
                   [exPAT, exCON] = se_GroupMatching(isPat,nPAT,nCON,Cov(:,ge),[Cov(:,ag) DVARS' FD' RMS']);
                    hasRS([exPAT exCON]) = 0;
                    hasRSfix([exPAT exCON]) = 0;
                else min(numel(nPAT),numel(nCON))>9;
                   [exPAT, exCON] = se_GroupMatching(isPat,nPAT,nCON,Cov(:,ge),Cov(:,ag));
                    hasRS([exPAT exCON]) = 0;
                    hasRSfix([exPAT exCON]) = 0;
                end
            else
                hasRS([nPAT nCON]) = 0;
                hasRSfix([nPAT nCON]) = 0;
            end
        end
        fprintf(1,'\n');
            
    end
end

if sum(hasRS)+sum(hasVBM)+sum(hasDTI)==0
    error('No subjects left after matching');
end


if who==1
    use = sum([hasRS; hasVBM; hasDTI])>0;
elseif who == 2
    use = (hasRS+hasVBM)==2;
else
    use = sum([hasRS; hasVBM; hasDTI])==3;
end

sub = sub(use);
SubDir = SubDir(use);
Cov = Cov(use,:);
hasRS  = hasRS(use);
hasRSfix  = hasRSfix(use);
hasVBM = hasVBM(use);
hasDTI = hasDTI(use);
isPat  = isPat(use);
scanner = scanner(use);

Cov(:,end+1) = tbv(use,1);  Covariates{end+1} = 'GMV';
Cov(:,end+1) = tbv(use,1)+tbv(use,2);  Covariates{end+1} = 'TBV';
Cov(:,end+1) = tbv(use,1)+tbv(use,2)+tbv(use,3);  Covariates{end+1} = 'ICV';
Cov(:,end+1) = GMmean(use);  Covariates{end+1} = 'GMmean';

Cov(:,end+1) = DVARS(use)';  Covariates{end+1} = 'DVARS';
Cov(:,end+1) = FD(use)';  Covariates{end+1} = 'FD';
Cov(:,end+1) = RMS(use)';  Covariates{end+1} = 'RMD';

QC = [DVARS(use)' FD(use)' RMS(use)'];
qc = {'DVARS', 'FD', 'RMS'};

try
    Qoff = find(isPat==1); Qon = [];
    for i=1:numel(Qoff);
        Qon = [Qon find(strcmpi(sub,strrep(sub(Qoff(i)),'_OFF','')))];
    end
    matchQ = [Qon; Qoff]';
catch
    matchQ = [];
end


save(fullfile(pwd,'Lookups',[spm_str_manip(gname,'rt') '_' suffix adfix '.mat']), ...
    'Covariates','sub','isPat','SubDir','Cov','hasRS','hasRSfix','hasDTI','hasVBM','matchQ');

fprintf(1,'\n%s\n',['Lookup created for ' spm_str_manip(gname,'rt')])

[~, ~, ~] = mkdir(fullfile(pwd,'Lookups','info')); 
fid = fopen(fullfile(pwd,'Lookups','info',[spm_str_manip(gname,'rt') '_' suffix adfix '.txt']),'w+');

fprintf(fid,'%s\n\n',['Sample information for Lookup-file: ' spm_str_manip(gname,'rt') '_' suffix adfix '.mat']);

for i=1:numel(sites)
    A = SubDir(scanner==sites(i));
    A = A{1};
    fprintf(fid,'%s\n',['Site ' int2str(i) ': ' A]);
end

hasThat = {hasRS hasVBM hasDTI};
whathas = {'RESTING-STATE data','VBM data','DTI data'};
ag = find(strcmpi(Covariates,'Age')==1);
ge = find(cellfun('length',strfind(lower(Covariates),'sex'))+ ...
    cellfun('length',strfind(lower(Covariates),'gender')));
ge = ge(1);

if any(isPat==1)
    
    for what = 1:3
        if mean(hasThat{what})>.2
            fprintf(fid,'\n%s\n\n%s\n%s\n\n%s\n\n','_________________________', ...
                whathas{what},['n = ' int2str(sum(hasThat{what}>0))],['PAT / CON: ' ...
                int2str(sum(isPat'==1 & hasThat{what}'>0)) ' / ' int2str(sum(isPat'==0 & hasThat{what}'>0))]);
            if numel(sites)>1          
                for i=1:numel(sites)
                    if sum(hasThat{what}'>0 & scanner==sites(i))>0
                        fprintf(fid,'%s\n',['Site ' int2str(i) ':  (PAT / CON)  ' ...
                            int2str(sum(isPat'==1 & hasThat{what}'>0 & scanner==sites(i))) ' / ' ...
                            int2str(sum(isPat'==0 & hasThat{what}'>0 & scanner==sites(i)))]);
                    end
                end
            end
            fprintf(fid,'\n\n%s\n\n','AGE comparison:');
            coh = {'Controls','Patients'};           
            for p = 0:1
                dat = Cov(isPat'==p & hasThat{what}'>0,ag);
                fprintf(fid,'%s%8.2f%8.2f%8.2f%8.2f\n',['All: ' coh{p+1} ' [mean, SD, median, IQR]'], ...
                    nanmean(dat),nanstd(dat),nanmedian(dat),iqr(dat));
            end
           
                [h p1] = ttest2(Cov(isPat'==0 & hasThat{what}'>0,ag),Cov(isPat'==1 & hasThat{what}'>0,ag));
                p2 = ranksum(Cov(isPat'==0 & hasThat{what}'>0,ag),Cov(isPat'==1 & hasThat{what}'>0,ag));
                fprintf(fid,'\t%s%10.4f%10.4f\n\n','ttest / ranksum : ' ,p1,p2);
           
            if numel(sites)>1
                for i=1:numel(sites)
                    if sum(hasThat{what}'>0 & scanner==sites(i))>0
                        for p = 0:1
                            dat = Cov(isPat'==p & hasThat{what}'>0 & scanner==sites(i),ag);
                            fprintf(fid,'%s%8.2f%8.2f%8.2f%8.2f\n',['Site ' int2str(i) ': '...
                                coh{p+1} '       [mean, SD, median, IQR]'],nanmean(dat),nanstd(dat),nanmedian(dat),iqr(dat));
                        end
                        
                            [h p1] = ttest2(Cov(isPat'==0 & hasThat{what}'>0 & scanner==sites(i),ag),...
                                Cov(isPat'==1 & hasThat{what}'>0 & scanner==sites(i),ag));
                            p2 = ranksum(Cov(isPat'==0 & hasThat{what}'>0 & scanner==sites(i),ag),...
                                Cov(isPat'==1 & hasThat{what}'>0 & scanner==sites(i),ag));
                            fprintf(fid,'\t%s%10.4f%10.4f\n','ttest / ranksum : ' ,p1,p2);
                        
                    end
                end
            end

            fprintf(fid,'\n\n%s\n\n','SEX comparison:');
                gc = unique(Cov(:,ge));
                if numel(gc)==2
                        for p = 0:1
                            dat = Cov(isPat'==p & hasThat{what}'>0,ge);
                            fprintf(fid,'%s%4.0f%4.0f\n',['All ' coh{p+1} ' [male, female]'],sum(dat==gc(1)), ...
                                sum(dat==gc(2)));
                        end
                    [table,chi2,p] = crosstab(isPat(hasThat{what}'>0)', Cov(hasThat{what}'>0,ge));
                    fprintf(fid,'\t%s%10.4f\n\n','chi2 : ' ,p);
                end
            if numel(sites)>1
                for i=1:numel(sites)
                    if sum(hasThat{what}'>0 & scanner==sites(i))>0
                        for p = 0:1
                            dat = Cov(isPat'==p & hasThat{what}'>0 & scanner==sites(i),ge);
                            fprintf(fid,'%s%4.0f%4.0f\n',['Site ' int2str(i) ': ' coh{p+1} '       [male, female]'], ...
                                sum(dat==gc(1)), sum(dat==gc(2)));
                        end
                        
                            [table,chi2,p] = crosstab(isPat(hasThat{what}'>0 & scanner==sites(i))', ...
                                Cov(hasThat{what}'>0 & scanner==sites(i),ge));
                        
                        fprintf(fid,'\t%s%10.4f\n','chi2 : ' ,p);
                    end
                end
            end

            
            if what==1
                for qci = 1:3
                    fprintf(fid,'\n\n%s\n\n',[qc{qci} ' comparison:']);
                    for p = 0:1
                        dat = QC(isPat'==p & hasThat{what}'>0,qci);
                        fprintf(fid,'%s%8.2f%8.2f%8.2f%8.2f\n',['All: ' coh{p+1} ...
                            ' [mean, SD, median, IQR]'],nanmean(dat),nanstd(dat),nanmedian(dat),iqr(dat));
                    end
                    [h p1] = ttest2(QC(isPat'==0 & hasThat{what}'>0,qci),QC(isPat'==1 & hasThat{what}'>0,qci));
                    p2 = ranksum(QC(isPat'==0 & hasThat{what}'>0,qci),QC(isPat'==1 & hasThat{what}'>0,qci));
                    fprintf(fid,'\t%s%10.4f%10.4f\n\n','ttest / ranksum : ' ,p1,p2);
                    if numel(sites)>1
                        for i=1:numel(sites)
                            if sum(hasThat{what}'>0 & scanner==sites(i))>0
                                for p = 0:1
                                    dat = QC(isPat'==p & hasThat{what}'>0 & scanner==sites(i),qci);
                                    fprintf(fid,'%s%8.2f%8.2f%8.2f%8.2f\n',['Site ' int2str(i) ': ' coh{p+1} ...
                                        '       [mean, SD, median, IQR]'],nanmean(dat),nanstd(dat),nanmedian(dat),iqr(dat));
                                end
                                
                                    [h p1] = ttest2(QC(isPat'==0 & hasThat{what}'>0 & scanner==sites(i),qci), ...
                                        QC(isPat'==1 & hasThat{what}'>0 & scanner==sites(i),qci));
                                    p2 = ranksum(QC(isPat'==0 & hasThat{what}'>0 & scanner==sites(i),qci), ...
                                        QC(isPat'==1 & hasThat{what}'>0 & scanner==sites(i),qci));
                                    fprintf(fid,'\t%s%10.4f%10.4f\n','ttest / ranksum : ' ,p1,p2);
                            
                            end
                        end
                    end
                end
                
                if (sum(hasRS) ~= sum(hasRSfix)) && sum(hasRSfix)>sum(hasRS)*0.7
                    fprintf(fid,'\n\n%s\n\n%s\n%s\n\n%s\n\n','------------------------',['Fixed RS sample is missing ' int2str(sum(hasRS) - sum(hasRSfix)) ...
                        ' subjects'],['n = ' int2str(sum(hasRSfix>0))],['PAT / CON: ' int2str(sum(isPat==1 & hasRSfix>0)) ...
                        ' / ' int2str(sum(isPat==0 & hasRSfix>0))],'AGE comparison:');
                    for p = 0:1
                        dat = Cov(isPat==p & hasRSfix>0,ag);
                        fprintf(fid,'%s%8.2f%8.2f%8.2f%8.2f\n',['All: ' coh{p+1} ' [mean, SD, median, IQR]'], ...
                            nanmean(dat),nanstd(dat),nanmedian(dat),iqr(dat));
                    end
                    
                        [h p1] = ttest2(Cov(isPat==0 & hasRSfix>0,ag),Cov(isPat==1 & hasRSfix>0,ag));
                        p2 = ranksum(Cov(isPat==0 & hasRSfix>0,ag),Cov(isPat==1 & hasRSfix>0,ag));
                        fprintf(fid,'\t%s%10.4f%10.4f\n\n','ttest / ranksum : ' ,p1,p2);
                  
                    if numel(sites)>1
                        for i=1:numel(sites)
                            if sum(hasRSfix'>0 & scanner==sites(i))>0
                                for p = 0:1
                                    dat = Cov(isPat'==p & hasRSfix'>0 & scanner==sites(i),ag);
                                    fprintf(fid,'%s%8.2f%8.2f%8.2f%8.2f\n',['Site ' int2str(i) ': ' coh{p+1} ...
                                        '       [mean, SD, median, IQR]'],nanmean(dat),nanstd(dat),nanmedian(dat),iqr(dat));
                                end
                                
                                    [h p1] = ttest2(Cov(isPat'==0 & hasRSfix'>0 & scanner==sites(i),ag), ...
                                        Cov(isPat'==1 & hasRSfix'>0 & scanner==sites(i),ag));
                                    p2 = ranksum(Cov(isPat'==0 & hasRSfix'>0 & scanner==sites(i),ag), ...
                                        Cov(isPat'==1 & hasRSfix'>0 & scanner==sites(i),ag));
                                    fprintf(fid,'\t%s%10.4f%10.4f\n','ttest / ranksum : ' ,p1,p2);
                                
                            end
                        end
                    fprintf(fid,'\n');
                    end

                    fprintf(fid,'\n%s\n','SEX comparison:');
                    if numel(gc)==2
                        for p = 0:1
                            dat = Cov(isPat'==p & hasRSfix'>0,ge);
                            fprintf(fid,'%s%4.0f%4.0f\n',['All ' coh{p+1} ' [male, female]'],sum(dat==gc(1)), ...
                                sum(dat==gc(2)));
                        end
                        [table,chi2,p] = crosstab(isPat(hasRSfix'>0)', Cov(hasRSfix'>0,ge));
                        fprintf(fid,'\t%s%10.4f\n','chi2 : ' ,p);
                    end
                    if numel(sites)>1
                    	for i=1:numel(sites)
                        	
                                if sum(hasRSfix'>0 & scanner==sites(i))>0
                                    for p = 0:1
                                        dat = Cov(isPat'==p & hasRSfix'>0 & scanner==sites(i),ge);
                                        fprintf(fid,'%s%4.0f%4.0f\n',['Site ' int2str(i) ': ' coh{p+1} '       [male, female]'], ...
                                            sum(dat==gc(1)), sum(dat==gc(2)));
                                    end
                                    [table,chi2,p] = crosstab(isPat(hasRSfix'>0 & scanner==sites(i))', ...
                                        Cov(hasRSfix'>0 & scanner==sites(i),ge));
                                    fprintf(fid,'\t%s%10.4f\n','chi2 : ' ,p);
                                end
                        	
                        end
                    end
                end
            end
        end
        fprintf(fid,'\n\n');
    end
    
else
    
    for what = 1:3
        if mean(hasThat{what})>.2
            fprintf(fid,'\n%s\n\n%s\n%s\n\n','_________________________',whathas{what},['n = ' int2str(sum(hasThat{what}>0))]);
            if numel(sites)>1          
                for i=1:numel(sites)
                    if sum(hasThat{what}'>0 & scanner==sites(i))>0
                        fprintf(fid,'%s\n',['Site ' int2str(i) ': ' int2str(sum(hasThat{what}'>0 & scanner==sites(i))) ...
                            ' subjects']);
                    end
                end
            end
            fprintf(fid,'\n\n%s\n\n','AGE comparison:');
                dat = Cov(hasThat{what}'>0,ag);
                fprintf(fid,'%s%8.2f%8.2f%8.2f%8.2f\n',['All subjects: [mean, SD, median, IQR]'], ...
                    nanmean(dat),nanstd(dat),nanmedian(dat),iqr(dat));
            if numel(sites)>1
                for i=1:numel(sites)
                    if sum(hasThat{what}'>0 & scanner==sites(i))>0
                            dat = Cov(hasThat{what}'>0 & scanner==sites(i),ag);
                            fprintf(fid,'%s%8.2f%8.2f%8.2f%8.2f\n',['Site ' int2str(i) ':       [mean, SD, median, IQR]'],...
                                nanmean(dat),nanstd(dat),nanmedian(dat),iqr(dat));
                     end
                end
            end

            fprintf(fid,'\n\n%s\n\n','SEX comparison:');
                gc = unique(Cov(:,ge));
                if numel(gc)==2

                            dat = Cov(hasThat{what}'>0,ge);
                            fprintf(fid,'%s%4.0f%4.0f\n',['All subjects: [male, female]'],sum(dat==gc(1)), ...
                                sum(dat==gc(2)));
                end
            if numel(sites)>1
                for i=1:numel(sites)
                 
                        if sum(hasThat{what}'>0 & scanner==sites(i))>0
                                dat = Cov(hasThat{what}'>0 & scanner==sites(i),ge);
                                fprintf(fid,'%s%4.0f%4.0f\n',['Site ' int2str(i) ':       [male, female]'], ...
                                    sum(dat==gc(1)), sum(dat==gc(2)));
                        end
                    
                end
            end

            
            if what==1
                for qci = 1:3
                    fprintf(fid,'\n\n%s\n\n',[qc{qci} ' comparison:']);
                    dat = QC(hasThat{what}'>0,qci);
                    fprintf(fid,'%s%8.2f%8.2f%8.2f%8.2f\n',['All subjects: [mean, SD, median, IQR]'],...
                        nanmean(dat),nanstd(dat),nanmedian(dat),iqr(dat));
                    if numel(sites)>1
                        for i=1:numel(sites)
                            if sum(hasThat{what}'>0 & scanner==sites(i))>0
                                dat = QC(hasThat{what}'>0 & scanner==sites(i),qci);
                                fprintf(fid,'%s%8.2f%8.2f%8.2f%8.2f\n',['Site ' int2str(i) ':       [mean, SD, median, IQR]'],...
                                    nanmean(dat),nanstd(dat),nanmedian(dat),iqr(dat));
                            end
                        end
                    end
                end
                
                if (sum(hasRS) ~= sum(hasRSfix)) && sum(hasRSfix)>sum(hasRS)*0.7
                    fprintf(fid,'\n\n%s\n\n%s\n%s\n\n%s\n\n','------------------------',['Fixed RS sample is missing ' int2str(sum(hasRS) - sum(hasRSfix)) ...
                        ' subjects'],['n = ' int2str(sum(hasRSfix>0))],'AGE comparison:');
                        dat = Cov(hasRSfix'>0,ag);
                        fprintf(fid,'%s%8.2f%8.2f%8.2f%8.2f\n',['All: subjects [mean, SD, median, IQR]'], ...
                            nanmean(dat),nanstd(dat),nanmedian(dat),iqr(dat));

                    if numel(sites)>1
                        for i=1:numel(sites)
                            if sum(hasRSfix'>0 & scanner==sites(i))>0
                                dat = Cov(hasRSfix'>0 & scanner==sites(i),ag);
                                fprintf(fid,'%s%8.2f%8.2f%8.2f%8.2f\n',['Site ' int2str(i) ':       [mean, SD, median, IQR]'],...
                                    nanmean(dat),nanstd(dat),nanmedian(dat),iqr(dat));
                             end
                        end
                    fprintf(fid,'\n');
                    end

                    fprintf(fid,'\n%s\n\n','SEX comparison:');
                    if numel(gc)==2
                            dat = Cov(hasRSfix'>0,ge);
                            fprintf(fid,'%s%4.0f%4.0f\n',['All subjects: [male, female]'],sum(dat==gc(1)), ...
                                sum(dat==gc(2)));
                     end
                    if numel(sites)>1
                    	for i=1:numel(sites)
                        	
                                if sum(hasRSfix'>0 & scanner==sites(i))>0
                                        dat = Cov(hasRSfix'>0 & scanner==sites(i),ge);
                                        fprintf(fid,'%s%4.0f%4.0f\n',['Site ' int2str(i) ':       [male, female]'], ...
                                            sum(dat==gc(1)), sum(dat==gc(2)));
                                end
                        	
                        end
                    end
                end
            end
        end
        fprintf(fid,'\n\n');
    end
    
end

fprintf(fid,'\n\n%s\n\n', '--- included subjects ---');

if exist('hasRSfix','var')
    fprintf(fid,'%s\t%s\t%s\t%s\t%s\n','subject','RS','RSfix','VBM','DWI');
    for i=1:numel(sub)
    fprintf(fid,'%s\t%s\t%s\t%s\t%s\n',sub{i},int2str(hasRS(i)),int2str(hasRSfix(i)), ...
        int2str(hasVBM(i)),int2str(hasDTI(i)));
    end
else
    fprintf(fid,'%s\t%s\t%s\t%s\n','subject','RS','VBM','DWI');
    for i=1:numel(sub)
        fprintf(fid,'%s\t%s\t%s\t%s\n',sub{i},int2str(hasRS(i)),int2str(hasVBM(i)),int2str(hasDTI(i)));
    end
end

fclose(fid);


if mean(hasRS)>.2
    load(fullfile(pwd,'matfiles','se_makeMean.mat')); cnt = 1;
    for s=1:numel(sub)
        if hasRS(s)
            matlabbatch{1}.spm.util.imcalc.input{cnt,1} = ...
                fullfile(datadir,'DATA',SubDir{s},sub{s},'RS',['c1w' sub{s} '_meanEPI.nii']);
            cnt = cnt+1;
        end
    end
    matlabbatch{1}.spm.util.imcalc.output = [spm_str_manip(gname,'rt') '_' suffix adfix '.nii'];
    matlabbatch{1}.spm.util.imcalc.outdir = {fullfile(pwd,'MaskenEtc','forRS')};
    spm_jobman('run',matlabbatch)

    load(fullfile(pwd,'matfiles','se_makeMean.mat')); cnt = 1;
    for s=1:numel(sub)
        if hasRS(s)
            matlabbatch{1}.spm.util.imcalc.input{cnt,1} = ...
                fullfile(datadir,'DATA',SubDir{s},sub{s},'RS',['w' sub{s} '_meanEPI.nii']);
            cnt = cnt+1;
        end
    end
    matlabbatch{1}.spm.util.imcalc.output = [spm_str_manip(gname,'rt') '_' suffix adfix '_meanEPI.nii'];
    matlabbatch{1}.spm.util.imcalc.outdir = {fullfile(pwd,'MaskenEtc','forRS')};
    spm_jobman('run',matlabbatch)
end

if mean(hasVBM)>.2
    load(fullfile(pwd,'matfiles','se_makeMean.mat')); cnt = 1;
    for s=1:numel(sub)
        if hasVBM(s)
            matlabbatch{1}.spm.util.imcalc.input{cnt,1} = ...
                fullfile(datadir,'DATA',SubDir{s},sub{s},'3D',['wmr' sub{s} '.nii']);
            cnt = cnt+1;
        end
    end
    matlabbatch{1}.spm.util.imcalc.output = [spm_str_manip(gname,'rt') '_' suffix adfix '_meanT1.nii'];
    matlabbatch{1}.spm.util.imcalc.outdir = {fullfile(pwd,'MaskenEtc','forVBM')};
    spm_jobman('run',matlabbatch)


    load(fullfile(pwd,'matfiles','se_makeMean.mat')); cnt = 1;
    for s=1:numel(sub)
        if hasVBM(s)
            matlabbatch{1}.spm.util.imcalc.input{cnt,1} = ...
                fullfile(datadir,'DATA',SubDir{s},sub{s},'3D',['m0wrp1' sub{s} '.nii']);
            cnt = cnt+1;
        end
    end
    matlabbatch{1}.spm.util.imcalc.output = [spm_str_manip(gname,'rt') '_' suffix adfix '_c1meanT1.nii'];
    matlabbatch{1}.spm.util.imcalc.outdir = {fullfile(pwd,'MaskenEtc','forVBM')};
    spm_jobman('run',matlabbatch)
end

fprintf(1,'\n%s\n',['All masks created for ' spm_str_manip(gname,'rt')]);


