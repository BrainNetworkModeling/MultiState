if ~(exist('VOIs','dir')==7)
  mkdir('VOIs');
end


VM = spm_vol(spm_select(Inf,'image','Select image files'));


mode = spm_input('Coordinates for','!+1','b','Max|Median|Peaks',[1 0 2],1,'Select mode');

u  = spm_input(['Hight threshold (0 for none)'],'+0','r',0,1);
k  = spm_input(['Extend threshold (0 for none)'],'+0','r',10,1);

if numel(VM)==1
    titlestr = spm_input('Project Name','+1','s',spm_str_manip(VM.fname,'rt'));
else
    titlestr = spm_input('Project Name','+1','s');
end


cnt = 1;
for bild=1:numel(VM)
    msk = spm_read_vols(VM(bild));
    ind = find(msk>u);
    [tX tY tZ] = ind2sub(VM(bild).dim,ind);
    pXYZ = [tX tY tZ]';
    Z = msk(ind);
    A    = spm_clusters(pXYZ);
    tXYZ = nan(3,max(A));
    
    
    fg = spm_figure('GetWin','Graphics'); spm_figure('Clear','Graphics'); spm_orthviews('Reset');
    spm_orthviews('Image', spm_vol(fullfile(spm('dir'),'canonical','single_subj_T1.nii')), [0.0 0.22 1 .8]);
    
    for cl=1:max(A)
        if sum(A==cl)>=k
            if mode == 1
                [voxel maxZ maxM dummy] = spm_max(Z(A == cl),pXYZ(:,A == cl));
                xtXYZ = maxM(:,find(maxZ==max(maxZ)));
            elseif mode == 0
                xtXYZ = round(median(pXYZ(:,A == cl)'))';
            else
                [voxel maxZ xtXYZ dummy] = spm_max(Z(A == cl),pXYZ(:,A == cl));
                [maxZ I] = sort(maxZ,'descend');
                xtXYZ    = xtXYZ(:,I(1:min(numel(I),3)));
            end
            
            
            AllmmXYZ = VM(bild).mat * [xtXYZ; ones(1,size(xtXYZ,2))];
            for rx = 1:size(xtXYZ,2)
                mmXYZ = AllmmXYZ(1:3,rx);
                spm_orthviews('addcolouredblobs',1,pXYZ(:,A==cl),Z(A==cl),VM(bild).mat,[1 0 0]);
                spm_orthviews('reposition',mmXYZ(1:3)')
                
                if numel(VM)>1
                    xVOI = spm_input(spm_str_manip(VM(bild).fname,'rt'),1,'s',spm_str_manip(VM(bild).fname,'rt'));
                else
                  if mode == 2
                      xVOI = spm_input(spm_str_manip(VM(bild).fname,'rt'),1,'s',num2str(maxZ(rx)));
                  else
                      xVOI = spm_input(spm_str_manip(VM(bild).fname,'rt'),1,'s',int2str(sum(A==cl)));
                    end
                end
                if ~strcmp(xVOI,'0')
                    mXYZ(:,cnt) = mmXYZ;
                    tXYZ(:,cnt) = xtXYZ(1:3,rx);
                    VOI{cnt}    = xVOI;
                    cnt = cnt+1;
                end
                spm_orthviews('rmblobs',1);
            end
        end
    end
    
end

mXYZ = mXYZ(1:3,:);



outp = 1;

if outp==1
    fid = fopen(fullfile(pwd,'VOIs',[titlestr '_VOIs.txt']),'wt');
    for i=1:numel(VOI)
        fprintf(fid,'%4.1f\t%3.1f\t%3.1f\t%s\n',mXYZ(:,i),VOI{i});
    end
    fclose(fid);
else
    save(fullfile(pwd,'VOIs',[titlestr '_VOIs.mat']),'VOI','tXYZ','VM','mXYZ');
end


OverGUI