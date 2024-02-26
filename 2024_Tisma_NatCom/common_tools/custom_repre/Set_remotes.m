function [inpath_bulk,inpath_group_test,outpath_group_test]= Set_remotes 
%this function sets local and remote paths depending where you run it

%it detects where I run it. local options are:
%laptop home
%laptop lab/office PC

inpath_bulk=[];
inpath_group_test=[];
outpath_group_test=[];


 %via  office PC:
if exist('N:\tnw\BN\CD\Shared','dir'),    
      inpath_bulk='N:\tnw\BN\CD\Shared\';
      inpath_group_test='M:\tnw\bn\cd\Shared\Jacob Kerssemakers\TESTdata_in\';
      outpath_group_test='M:\tnw\bn\cd\Shared\Jacob Kerssemakers\TESTdata_out\';
end 
 
%via laptop & VPN:
 if exist('V:\tnw\bn\cd\Shared\','dir')   
      inpath_bulk='X:\tnw\BN\CD\Shared\';
      inpath_group_test='V:\tnw\bn\cd\Shared\Jacob Kerssemakers\TESTdata_in\';      
      outpath_group_test='V:\tnw\bn\cd\Shared\Jacob Kerssemakers\TESTdata_out\';      
 end   
 
 
 