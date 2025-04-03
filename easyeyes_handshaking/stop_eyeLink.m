Eyelink('stoprecording');
mytable = struct2table(dataEyeLink);
save(sprintf('%s/%s.mat',oo.dataFolder,oo.dataFilename),'oo','mytable');

Eyelink('command','clear_screen');
Eyelink('command', 'record_status_message ''END''');

statRecFile = Eyelink('ReceiveFile',const.edffilename,const.edffilename);

if statRecFile ~= 0
    fprintf(1,'\n\tEyelink EDF file correctly transferred');
else
    fprintf(1,'\n\Error in Eyelink EDF file transfer');
    statRecFile2 = Eyelink('ReceiveFile',const.edffilename,const.edffilename);
    if statRecFile2 == 0
        fprintf(1,'\n\tEyelink EDF file is now correctly transferred');
    else
        fprintf(1,'\n\n\t!!!!! Error in Eyelink EDF file transfer !!!!!');
        my_sound(9,aud);
    end
end
Eyelink('CloseFile');
WaitSecs(2.0);
Eyelink('Shutdown');
WaitSecs(2.0);
