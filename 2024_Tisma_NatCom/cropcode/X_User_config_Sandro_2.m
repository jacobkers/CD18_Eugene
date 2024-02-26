 function pathinfo=X_User_config_Sandro_2(expno);
        pathinfo.maxcellsperframe=10E6;  %limited file run for test purposes; set to 1E9 to run all data
        pathinfo.autorun=1;
 
        pathinfo.dircode = 'C:\Users\aleksandrejapa\Desktop\CellCropCode'; % Directory of the codes
        pathinfo.dirdipimage = 'C:\Program Files\DIPimage 2.9\common\dipimage'; % Directory of dipimage         
        switch expno
            case 1
               pathinfo.experimentpath = 'C:\Data\Hi-C samples\2018.08.20-2830 MG1655\30C\2830+30C-OD_02.001\';
               pathinfo.centerplane = [3,3]; % Needs review in fiji!!
               pathinfo.limitzplanes=1;
               pathinfo.txycz_template='2830+30C-OD_02.001xy1c1z1.tif';
        end