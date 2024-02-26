function A025_content_plots(init,expi,usr)
%JWJK_A:-------------------------------------------------------------------
%Summary: %This function analyzes spots positions associated with 
%DNA plectonemes and condensin
%Approach: the positions of condensin and plectonemes are
%related to each other; peaks are counted.
%References: CD lab, project Eugene Kim, written by Jacob Kers, 2019
%:JWJK_A-------------------------------------------------------------------
close all;

%% 1) Set common paths; use standardized naming

AllExp=init.AllExp;
inpath=strcat(init.datapathout, 'matlabresults\',init.expname,'\');
outpath=strcat(init.datapathout, 'matlabresults\',init.expname,'A025_contentplots\');
if ~isdir(outpath), mkdir(outpath); end;
LE=length(AllExp);  %for all experiments
for roi=1:LE  
    Exp=strcat(init.roidirname,num2str(AllExp(roi)));
    roino=AllExp(roi);
    switch usr
    case 'Jacob',  expinfo=A002_JK_Condensin_with_plectonemes_expinfo(expi,AllExp(roi));
    case 'Eugene', expinfo=A002_EK_Condensin_with_plectonemes_expinfo(expi,AllExp(roi));
    end    
    
    LoadName=char(strcat(inpath, 'EKMcp_A020_',Exp));             
    disp(strcat('A025_Contentplots: Exps to work through:',num2str(LE-roi)));
    load(strcat(LoadName, '_allresults.mat'));     
    [kymo_duration,kymo_width]=size(kymo_DNA); 
    
    SaveName=[outpath,'EKMcp_A025_roi', num2str(roi)];
    content_of_plectonemes(info_DNA,kymo_DNA,SaveName);
    dum=1;
end 



  



        