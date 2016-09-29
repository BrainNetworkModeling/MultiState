function ar_SelectSubjects_nogui(path,list_file)

% Creates a list of subjects for preprocessing.
% 
% path 			 - folder of the cohort to process
% list_file      - name of text file containing a list of subjects to process



%--------------- readout of path-to-data ----------------------------
fid = fopen('_path2data','r'); datadir = fgetl(fid); fclose(fid);


warning off all;

if nargin == 0
    path = input('Enter folder name of the cohort to preprocess (''sample''): ');
    list_file = input('Enter name of subject list file (''sample.txt''): ');
end

names = dir(fullfile(datadir,'DATA',path));
names = names([names.isdir]); names = names(3:end);


fid = fopen(fullfile(pwd,list_file),'w+');

for i = 1 : size(names,1)
    fprintf(fid, '%s\n', fullfile(datadir,'DATA',path,names(i).name));
end

fclose(fid);

fprintf('List written to %s\n', fullfile(pwd,list_file));

end