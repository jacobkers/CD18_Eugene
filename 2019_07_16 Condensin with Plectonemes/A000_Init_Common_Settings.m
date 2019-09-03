function init=A000_Init_Common_Settings
%Settings applied to all files

%local paths
JKpth='C:\Users\jkerssemakers\'
init.datapathin=[JKpth 'CD_Data_in\2018_Eugene\'];                      %loading
init.datapathout=[JKpth 'Dropbox\CD_Data_out\2018_Eugene\'];            %saving
addpath(genpath([JKpth 'Dropbox\CD_recent\BN_CD18_Eugene\Matlab_tools_CD18EK'])); %tools:

init.psf_est=2.7; 
%estimated point spread function (used for estimated peak content and localization)

%kymograph settings
init.t_smooth=10;           %for kymographs, in frames 
init.x_smooth=4;            %for kymographs, in pixels 
init.tresholdsigmas=2;      %number of sigmas beyond background noise,..
                            %..used for getting fluorescence counts and levels


                          