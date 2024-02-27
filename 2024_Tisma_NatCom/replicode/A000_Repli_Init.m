function initval=A000_Repli_Init(batchrunindex,usr);
%JWJK_A:-------------------------------------------------------------------
%Initialization

%Description: set main paths of data to be analyzed. Settings per
%experiment and initialization 

%Reference: CD lab, project Sandro, written by Jacob Kers 2018-20
%:JWJK_A-------------------------------------------------------------------


switch usr
    case 'Sandro',initval=A_User_config_Sandro(batchrunindex);
    case 'Jacob' ,initval=A_User_config_Jacob(batchrunindex);     
end


initval.codepth=pwd;
addpath(initval.codepth);
addpath(strcat(initval.codepth,initval.DirSep,'tools_repli',initval.DirSep));
cd ..;
addpath(genpath(strcat(pwd,initval.DirSep,'common_tools',initval.DirSep))); 
cd(initval.codepth);


%check the crop code naming format. Only do this if the fields were not
%specified  before in config settings. Older runs do not store these
%fields!
if ~isfield(initval,'channelorder'), 
    [channelorder, ~]=Check_channel_order(initval); 
    initval.channelorder=channelorder;
end
if ~isfield(initval,'numberofchannels')
    [~, numberofchannels]=Check_channel_order(initval); 
    initval.numberofchannels=numberofchannels;
end


%% overwrite detailed info for experiment via Excel Table
if ~strcmp(initval.pth_excelfile, 'none'), initval=fetch_excel_repli(initval);  end  
initval.Cell_Labels=Select_data(initval); 
initval.blur_ter=7; 
initval=orderfields(initval);                       
                        
 function [channelorder,numberofchannels]=Check_channel_order(initval,epth);   
     %This fuction checks whci channel order was used    
     % OLD               NEW
    		% cellc1: cfp-ter   ph/bf
    		% cellc2: yfp-chro  chro (same)
    		% cellc3: rfp-ori   ori (same)
    		% cellc4: ph/bf     ter (cfp)
    		% cellc5: extra     mukbef (optional)
		%Thus,  c1 and c4 are effectively swapped!
		%note:	pathinfo.newchannelnaming is saved in X020: celldata.mat
    %'OLD'
    channelorder=1;
    numberofchannels=4;
    pth=strcat(initval.pth_crop, 'X020_cellcoordinate\');
    if isdir(pth)
    info=load([pth,'celldata.mat']);
    if isfield(info, 'pathinfo')
        %NEW'
        if isfield(info. pathinfo, 'newchannelnaming')
            channelorder=2;  
            numberofchannels=info.pathinfo.numberofcolours;
        end
    end
    end
    
 function initval=fetch_excel_repli(initval);
        [numdat,textdat]= xlsread(initval.pth_excelfile);
         Headers=textdat(1,:);
         col_Exps=find(strcmp(Headers,'ExpLabel')); %identify type column
         Exps=textdat(2:end,col_Exps);
         Exp_Row=find(strcmp(Exps,initval.expi)); 

         %get text fields from excel
         initval.searchlabel=char(textdat(Exp_Row+1,(find(strcmp(Headers,'searchlabel')))));
         initval.cellmaskname=char(textdat(Exp_Row+1,(find(strcmp(Headers,'MaskName')))));

         %get numeric parameters
         initval.analyze_ori_ter=(numdat(Exp_Row,(find(strcmp(Headers,'ori_ter')))-1)); 
         initval.analyze_replisome=(numdat(Exp_Row,(find(strcmp(Headers,'replisome')))-1)); 
         initval.Psf_est=(numdat(Exp_Row,(find(strcmp(Headers,'Pointspread')))-1)); 
         initval.cluster_sep_sigs=(numdat(Exp_Row,(find(strcmp(Headers,'SeparateSigmas')))-1)); 
         initval.nmperpixel=(numdat(Exp_Row,(find(strcmp(Headers,'NmPerPixel')))-1)); 
