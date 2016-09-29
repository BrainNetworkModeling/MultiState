
Di = spm_select(Inf,'dir','Select subject directories');

[list_file,path] = uiputfile('*.*','Select output file');

format = '%s\n';
fid = fopen(fullfile(path,list_file),'w+');

for i = 1 : size(Di,1)
    fprintf(fid, '%s\n', Di(i,:));
end

fclose(fid);

fprintf('List written to %s\n', fullfile(path,list_file));