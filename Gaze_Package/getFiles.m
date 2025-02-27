function [cursorFiles,mainFiles, eyelinkFiles] = getFiles(mydir)
% loads all files in a folder and categorize them into cursor files, main
% output (psychoJS files) from EasyEyes, and eyelink files.
% throws an error message when there are incomplete groups (e.g., only
% cursor file, but no main or eyelink files for a participant)

    d = dir(sprintf('%s/*.csv',mydir));
    
    files = {d.name};
    
    for f = 1 :length(files)
        
        cursor = dir(sprintf('%s/*%s*_cursor.csv',mydir,files{f}(1:3)));
        cursorFiles{f} = cursor.name;
    
        mainFile = dir(sprintf('%s/*%s*_main.csv',mydir,files{f}(1:3)));
        mainFiles{f} = mainFile.name;
    
        eyelinkFile = dir(sprintf('%s/*%s*_matlabOutput.csv',mydir,files{f}(1:3)));
        eyelinkFiles{f} = eyelinkFile.name;
    end
    
    cursorFiles = unique(cursorFiles);
    mainFiles = unique(mainFiles);
    eyelinkFiles = unique(eyelinkFiles);
    
    assert(length(cursorFiles) == length(mainFiles));
    assert(length(eyelinkFiles) == length(cursorFiles));

end