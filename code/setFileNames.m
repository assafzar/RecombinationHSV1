% dname = '/project/bioinformatics/Danuser_lab/shared/assaf/OrenKobilerTAU/201702/'
% addpath(genpath('/project/bioinformatics/Danuser_lab/shared/assaf/OrenKobilerTAU/code'));
function [] = setFileNames(dname)

filenames = dir(dname);
nfiles = length(filenames);

for i = 1 : nfiles
    filename = filenames(i).name;
    [pathstr, name, ext] = fileparts(filename);
    
    %% tif file
    if (strcmp(ext,'.tif'))
        if ~isempty(regexp(name,'[0-9]C[0-9]','match')) || ~isempty(regexp(name,'[0-9]C[0-9][0-9]','match'))
            continue;
        end
        curFname = [dname filesep name ext]; 
        if strcmp(name(1:3),'31 ')
            name = name(4:end);
        end
        name0 = strrep(name,'(','');
        name0 = strrep(name0,')','');
        name0a = strrep(name0,'vero 1','');
        name1 = strrep(name0a,' ','');
        name2 = strrep(name1,'.tif','');
        name3 = strrep(name2,'-','');
        name4 = strrep(name3,'.tif','');% in case some filenames are without the .tif
        name5 = strrep(name4,'=','');
        name6 = strrep(name5,'Vero31','');
        name7 = strrep(name6,'Vero','');        
        name7 = strrep(name7,'vero','');
        name8 = strrep(name7,'32X35','');
        name9 = strrep(name8,'25X35','');
        outStr = regexp(name9,'slice [0-9]','match');
        if ~isempty(outStr)
            name10 = strrep(name9,outStr{1},'');
        else
            outStr = regexp(name9,'slice[0-9]','match');
            if ~isempty(outStr)
                name10 = strrep(name9,outStr{1},'');
            else
                name10 = name9;
            end
        end
        if strcmp(name10,'C0') || strcmp(name10,'C1') || strcmp(name10,'C2')
            name10 = ['0C' name10(2)];
        end
        
        newFname = [dname filesep name10 ext]; 
        
        
        movefile(curFname,newFname);
    end
end
end
