function initval=A000_Lociicode_Get_PathsandExperiments(expstring);
%JWJK_B:-------------------------------------------------------------------
%
%Summary: Initialization of settings and paths of use for locii code
%
%References: code writen by JWJK,2017
%
%:JWJK_B-------------------------------------------------------------------
%bulk'and remote access

initval.pix2nm=60;
initval.shortset=10E6; %for testing, otherwise set to 10E6
initval.cropping_edge=20;  %for ROIs


if nargin<1, expstring=2.1;end
initval.runmodus='autorun'  ; %'user_screen'
%initval.runmodus='user_screen'  ; %



%% path settings, 
curpth=pwd;
cd .. %repo level
repopath=pwd;
addpath([repopath, '/common_tools']);
cd(curpth);


initval.codepth=[pwd '\'];

[inpath_bulk,inpath_group_test,outpath_group_test]= Set_remotes;
  

switch floor(expstring)
   case 1,initval.expname='2017_08_31_2179-rod-dynamics';
         initval.supermainpath='C:\jkerssemakers\My Documents\BN CD Data\2016_Sandro\';
         expth='2017_08_31_2179-rod-dynamics\';         
         switch expstring
             case 1.1, initval.subdir='2179-rod-1sec-dynamics001';
         end
         initval.imageindir='DriftCorrected\';
    case 2
        initval.expname='Locii Revisit';   
        initval.supermainpath=[inpath_bulk 'Timo Muller\']; 
        initval.supermainpath=['O:\CD_Data_in\2016_Sandro\' 'Timo Muller\']; %blue O:
         expth='2019.12.5 - 2179 room temp\';        
         initval.SaveLabels=[{'bf'},{'color1'},{'color2'},{'color3'} {'lab'}];
         %five elements, use 'X' for no channel
        
         switch expstring
             case 2.1, initval.subdir='2179_RT_001\';
             case 2.2, initval.subdir='2179_RT_002\';
             case 2.3, initval.subdir='2179_RT_003\';
             case 2.4, initval.subdir='2179_RT_004\'; 
             case 2.5, initval.subdir='2179_RT_005\';         
         end

         initval.imageindir='DriftCorrected\'; 
end


initval.maininpath=strcat(initval.supermainpath,expth);
initval.imagepath=strcat(initval.supermainpath,expth,initval.subdir);
initval.matfilesinoutdir=strcat(initval.imagepath,'matlabresults\');
if ~isfolder(initval.matfilesinoutdir), mkdir(initval.matfilesinoutdir); end
