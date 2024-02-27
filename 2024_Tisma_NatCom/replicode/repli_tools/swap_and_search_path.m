function outpth=swap_and_search_path(inpth, sourcetype)
%work with lokal dropboxes on different
%machines and not change too much old code.
%paths to check:

%1) standard work directories office PC or laptop
%one contains 'CD_Data_in'; 
%for data-in paths, should be next to Dropbox
%all others contain 'Dropbox'; 

%2) examples of (data) remote paths

%M:\tnw\bn\cd\Shared    %group drive office
%N:\tnw\BN\CD\Shared    %bulk drive office PC


%V:\tnw\bn\cd\Shared    %group drive remote VPN in lockdown times
%X:\tnw\BN\Shared\      %bulk drive remote VPN in lockdown times


%Jacob Kers, 2019
clc
if nargin<2
    sourcetype = 'project'; %'data_in' %'data_out'
    inpth='D:\flappaflapa\CD_Data_in\2016_Sandro'; 
    inpth='D:\trala\Dropbox\BN_recent\BN_All_14_Matlab_Code\JK_CleanTools';
    inpth='N:\tnw\BN\CD\Shared\Sandro\MukB-MatP paperN:\tnw\BN\CD\Shared\Sandro\MukB-MatP paper';
    inpth='X:\tnw\BN\Shared\Sandro\MukB-MatP paperN:\tnw\BN\CD\Shared\Sandro\MukB-MatP paper';
    inpth='2179-initiations_A1-1\';
    display(inpth);
end

%% 1) check what local path
    %outpth=Patch_Path(inpth,'Dropbox','dropbox_project');    
    outpth=Patch_Path(inpth,'CD_Data_in','data_in_search');
    outpth=Patch_Path(outpth,'CD_Data_out','data_output_search'); 

   dum=1;
    
    function outpth=Patch_Path(inpth,commondir,sourcetype);     
    outpth=inpth;                %            
    %% 1) first, adjust remote mapping
    inpth=adjust_remotes(inpth);
    
        
    %% next, find a local path containing common part
    switch sourcetype
        case 'dropbox_project'  %local path with dropbox
            startpth=char(pwd);            
        case 'data_in_search' %search work in order of preference

            %local laptop:
            %dir_list=find_data_dir('D:\jkerssemakers\CD_Data_in',inpth,1);
            %if isempty(dir_list)
                dir_list=find_data_dir(adjust_remotes('N:\tnw\BN\CD\Shared\Jacob Kerssemakers\CD_Data_in\'),inpth,1);
            %end
            %local O:
            %Kdrive testdata
            %bulk
            if exist('D:\jkerssemakers\')  %office PC
                startpth='D:\jkerssemakers\Dropbox\CD_Data_out';   
            end                            
    end

            k1 = strfind(startpth,char(commondir));
            %the path above Dropbox:
            if ~isempty(k1), upperpath=startpth(1:k1-1); end 
    
    if ~isempty(k1)
        k2 = strfind(inpth,char(commondir));
        if ~isempty(k2)
            lowerpath=inpth(k2:end);
            outpth=[upperpath lowerpath];
        else
            outpth=inpth;                %no idea
        end  
    end
    

    



