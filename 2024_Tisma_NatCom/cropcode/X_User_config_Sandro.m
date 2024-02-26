function pathinfo=X_User_config_Sandro(expno,override);
        
        pathinfo.maxcellsperframe=10E6;  %limited file run for test purposes; set to 1E9 to run all data       
        pathinfo.cropedge=5;        % cropping space in all directions
        pathinfo.autorun=0;

        pathinfo.dircode = pwd; % Directory of the codes
        pathinfo.dirdipimage = 'C:\Program Files\DIPimage 2.9\common\dipimage'; % Directory of dipimage
        pathinfo.mainpath='O:\My Nguyen\2820\2018-06-14-2820-DAPI\';
        
        if override
            pathinfo.mainpath='M:\tnw\bn\cd\Shared\Jacob Kerssemakers\TESTdata_in\Sandro\';
        end
        
        if nargin<1, expno=1;end
          switch expno
              case 0   %TEST     
                    pathinfo.experimentpath=[pathinfo.mainpath,'CropCodeTestData\'];
                    pathinfo.centerplane = [6 6 6]; 
                    pathinfo.limitzplanes=0;
                    % Center plane of z-scan in different channels
                     %3 entries corresponds to [brightfield, yfp, rfp]
                    pathinfo.txycz_template='2179-initiations_A1-1t01xy1c1z1.tif';
                case 1   
                    pathinfo.experimentpath=[pathinfo.mainpath '\2820\2018-06-14-2820-DAPI\2018-06-14-2820-DAPI\'];
                    pathinfo.centerplane = [8,7,7,8] ;
                             % Center plane of z-scan in different channels
                             %3 entries correspond to [brightfield, yfp, rfp]
                             %4 entries correspond to [brightfield, yfp, rfp cfp]
                             pathinfo.limitzplanes=1;
                    pathinfo.txycz_template='2018-06-14-2820-DAPI-3xy1c1z01.tif'; 
                    
          end