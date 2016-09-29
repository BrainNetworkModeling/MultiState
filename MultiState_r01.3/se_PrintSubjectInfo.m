clear
gname = spm_select(Inf,'any','Select lookup speadsheets','',fullfile(pwd,'Lookups'),'.mat');

load(gname);



clc

si = find(strcmpi(Covariates,'Site')==1);
sites = unique(Cov(:,si));

ag = find(strcmpi(Covariates,'Age')==1);
ge = find(cellfun('length',strfind(lower(Covariates),'sex'))+cellfun('length',strfind(lower(Covariates),'gender'))); ge = ge(1);

hasThat = {hasRS hasVBM hasDTI};
whathas = {'forRS','forVBM','forDTI'};

QC = [Cov(:,end-2)  Cov(:,end-1)  Cov(:,end) ];
qc = {'DVARS', 'FD', 'RMS'};

if any(isPat==1)
    
    fprintf(1,'\n\n\n\n\n\n')
    
    hasThat = {hasRS hasVBM hasDTI};
    whathas = {'forRS','forVBM','forDTI'};
    ag = find(strcmpi(Covariates,'Age')==1);
    ge = find(cellfun('length',strfind(lower(Covariates),'sex'))+cellfun('length',strfind(lower(Covariates),'gender'))); ge = ge(1);
    
    for what = 1:3
        if mean(hasThat{what})>.2
            fprintf(1,'\n\n\n%s\n\n',whathas{what});
            coh = {'Controls','Patients'};
            if any(strcmp(Covariates,'Site'))
                si = find(strcmpi(Covariates,'Site')==1);
                sites = unique(Cov(:,si));
            else
                Cov(:,end+1) = 1;
                si = size(Cov,2);
                sites = unique(Cov(:,si));
            end
            for i=1:numel(sites)
                if sum(hasThat{what}'>0 & Cov(:,si)==sites(i))>0
                    for p = 0:1
                        dat = Cov(isPat'==p & hasThat{what}'>0 & Cov(:,si)==sites(i),ag);
                        fprintf(1,'%s%8.2f%8.2f%8.2f%8.2f\n',['Site ' int2str(i) ': ' coh{p+1} ' [mean, SD, median, IQR]'],nanmean(dat),nanstd(dat),nanmedian(dat),iqr(dat))
                    end
                    try
                        [h p1] = ttest2(Cov(isPat'==0 & hasThat{what}'>0 & Cov(:,si)==sites(i),ag),Cov(isPat'==1 & hasThat{what}'>0 & Cov(:,si)==sites(i),ag));
                        p2 = ranksum(Cov(isPat'==0 & hasThat{what}'>0 & Cov(:,si)==sites(i),ag),Cov(isPat'==1 & hasThat{what}'>0 & Cov(:,si)==sites(i),ag));
                        fprintf(1,'\t%s%10.4f%10.4f\n','ttest / ranksum : ' ,p1,p2)
                    end
                end
            end
            fprintf(1,'\n');
            for p = 0:1
                dat = Cov(isPat'==p & hasThat{what}'>0,ag);
                fprintf(1,'%s%8.2f%8.2f%8.2f%8.2f\n',['All: ' coh{p+1} ' [mean, SD, median, IQR]'],nanmean(dat),nanstd(dat),nanmedian(dat),iqr(dat))
            end
            [h p1] = ttest2(Cov(isPat'==0 & hasThat{what}'>0,ag),Cov(isPat'==1 & hasThat{what}'>0,ag));
            p2 = ranksum(Cov(isPat'==0 & hasThat{what}'>0,ag),Cov(isPat'==1 & hasThat{what}'>0,ag));
            fprintf(1,'\t%s%10.4f%10.4f\n\n\n','ttest / ranksum : ' ,p1,p2)
            
            
            fprintf(1,'\n');
            for i=1:numel(sites)
                if sum(hasThat{what}'>0 & Cov(:,si)==sites(i))>0
                    gc = unique(Cov(:,ge));
                    for p = 0:1
                        dat = Cov(isPat'==p & hasThat{what}'>0 & Cov(:,si)==sites(i),ge);
                        fprintf(1,'%s%4.0f%4.0f\n',['Site ' int2str(i) ': ' coh{p+1} ],sum(dat==gc(1)), sum(dat==gc(2)))
                    end
                    [table,chi2,p] = crosstab(isPat(hasThat{what}'>0 & Cov(:,si)==sites(i))', Cov(hasThat{what}'>0 & Cov(:,si)==sites(i),ge));
                    fprintf(1,'\t%s%10.4f\n','chi2 : ' ,p)
                end
            end
            fprintf(1,'\n');
            for p = 0:1
                dat = Cov(isPat'==p & hasThat{what}'>0,ge);
                fprintf(1,'%s%4.0f%4.0f\n',['All ' coh{p+1} ],sum(dat==gc(1)), sum(dat==gc(2)))
            end
            [table,chi2,p] = crosstab(isPat(hasThat{what}'>0)', Cov(hasThat{what}'>0,ge));
            fprintf(1,'\t%s%10.4f\n','chi2 : ' ,p)
            
            
            fprintf(1,'\n');
            for i=1:numel(sites)
                if sum(hasThat{what}'>0 & Cov(:,si)==sites(i))>0
                    fprintf(1,'%s\n',['Site ' int2str(i) ':  (PAT / CON)  ' ...
                        int2str(sum(isPat'==1 & hasThat{what}'>0 & Cov(:,si)==sites(i))) ' / ' ...
                        int2str(sum(isPat'==0 & hasThat{what}'>0 & Cov(:,si)==sites(i)))])
                end
            end
            fprintf(1,'\n');
            fprintf(1,'\n%s\n\n',['All:  (PAT / CON)  ' ...
                int2str(sum(isPat'==1 & hasThat{what}'>0)) ' / ' ...
                int2str(sum(isPat'==0 & hasThat{what}'>0))])
            
            
            
            
            
            if what==1
                for qci = 1:3
                    fprintf(1,'\n%s\n',qc{qci})
                    for i=1:numel(sites)
                        if sum(hasThat{what}'>0 & Cov(:,si)==sites(i))>0
                            for p = 0:1
                                dat = QC(isPat'==p & hasThat{what}'>0 & Cov(:,si)==sites(i),qci);
                                fprintf(1,'%s%8.2f%8.2f%8.2f%8.2f\n',['Site ' int2str(i) ': ' coh{p+1} ' [mean, SD, median, IQR]'],nanmean(dat),nanstd(dat),nanmedian(dat),iqr(dat))
                            end
                            try
                                [h p1] = ttest2(QC(isPat'==0 & hasThat{what}'>0 & Cov(:,si)==sites(i),qci),QC(isPat'==1 & hasThat{what}'>0 & Cov(:,si)==sites(i),qci));
                                p2 = ranksum(QC(isPat'==0 & hasThat{what}'>0 & Cov(:,si)==sites(i),qci),QC(isPat'==1 & hasThat{what}'>0 & Cov(:,si)==sites(i),qci));
                                fprintf(1,'\t%s%10.4f%10.4f\n','ttest / ranksum : ' ,p1,p2)
                            end
                        end
                    end
                    for p = 0:1
                        dat = QC(isPat'==p & hasThat{what}'>0,qci);
                        fprintf(1,'%s%8.2f%8.2f%8.2f%8.2f\n',['All: ' coh{p+1} ' [mean, SD, median, IQR]'],nanmean(dat),nanstd(dat),nanmedian(dat),iqr(dat))
                    end
                    [h p1] = ttest2(QC(isPat'==0 & hasThat{what}'>0,qci),QC(isPat'==1 & hasThat{what}'>0,qci));
                    p2 = ranksum(QC(isPat'==0 & hasThat{what}'>0,qci),QC(isPat'==1 & hasThat{what}'>0,qci));
                    fprintf(1,'\t%s%10.4f%10.4f\n\n','ttest / ranksum : ' ,p1,p2)
                    
                    
                end
            end
        end
        fprintf(1,'\n\n')
    end
    
    
    for i=1:numel(sites)
        A = SubDir(Cov(:,si)==sites(i)); A = A{1};
        fprintf(1,'%s\n',['Site ' int2str(i) ': ' A])
    end
    
    
end

