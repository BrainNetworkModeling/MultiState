clear, clc; warning off all; spm('Defaults','fmri')


Di = spm_select(Inf,'dir','Select subject directories');

fid = spm_input('Output File','+1','s',['Movement_' date]);
fid = fopen(fullfile(pwd,[fid '.txt']),'w+');

fprintf(fid,'%s\t%s\t%s\t%s\n','Subject','RMS','FD','DVARS');

for i=1:size(Di,1)
    subjects{i} = Di(i,:);
    nx = Di(i,:);
    if nx(end) == filesep; nx(end) = []; end
    us = strfind(nx,filesep);
    nx = Di(i,us(end)+1:end);
    if nx(end) == filesep; nx(end) = []; end
    sub{i} = nx;
end

RMS   = nan(1,numel(sub));
FD    = nan(1,numel(sub));
DVARS = nan(1,numel(sub));

for i=1:numel(sub)
    if numel(dir(fullfile(subjects{i},'RS','ICA',[sub{i} '.nii.gz'])))
        info = load(fullfile(subjects{i},'RS',['Counfounds_' sub{i} '.mat']));
        rp = info.rp;
        
        drp = [zeros(1,6); diff(rp)];
        move = sqrt(sum(drp(:,1:3).^2,2));
        
        euler = zeros(size(rp,1),1);
        for ii=1:size(rp,1)
            euler(ii) = 50*acos((cos(drp(ii,4))*cos(drp(ii,5)) + cos(drp(ii,4))*cos(drp(ii,6)) + cos(drp(ii,5))*cos(drp(ii,6)) + sin(drp(ii,4))*sin(drp(ii,6))*sin(drp(ii,5)) - 1)/2);
        end
        
        RMS(i)  = sqrt(nanmean([move+euler].^2));
        tmp   = sum(abs([drp(:,1:3) drp(:,4:6)*50]),2);
        FD(i)   = sqrt(nanmean(tmp.^2));
        
        dat = spm_read_vols(spm_vol(fullfile(subjects{i},'RS',['MNI152GM_' sub{i} '.nii'])));
        dat = (diff(dat)./dat(1:end-1,:))*100;
        dat(abs(dat)>10) = NaN;
        DVARS(i) = sqrt(nanmean(nanmean(dat.^2,2)));
        
        fprintf(fid,'%s\t%6.5f\t%6.5f\t%6.5f\n',sub{i},RMS(i),FD(i),DVARS(i));
        fprintf(1,'%s\n',[int2str(i) ' / ' int2str(numel(sub))]);
    end
end