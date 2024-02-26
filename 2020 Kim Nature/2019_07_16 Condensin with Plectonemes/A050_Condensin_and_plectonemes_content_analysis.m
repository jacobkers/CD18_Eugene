function A050_Condensin_and_plectonemes_content_analysis(init,expi,usr)
%JWJK_A:-------------------------------------------------------------------
%Summary: %This function analyzes spots positions associated with 
%DNA plectonemes and condensin
%Approach: the positions of condensin and plectonemes are
%related to each other; peaks are counted.
%References: CD lab, project Eugene Kim, written by Jacob Kers, 2019
%:JWJK_A-------------------------------------------------------------------
close all;
actions.collectresults=1;
%% 1) Set common paths; use standardized naming

AllExp=init.AllExp;
inpath=strcat(init.datapathout, 'matlabresults\',init.expname,'\');
outpath=strcat(init.datapathout, 'matlabresults\',init.expname,'A050_contentplots\');
if ~isdir(outpath), mkdir(outpath); end;
if actions.collectresults
LE=length(AllExp);  %for all experiments
for roi=1:LE  
    Exp=strcat(init.roidirname,num2str(AllExp(roi)));
    roino=AllExp(roi);
    switch usr
    case 'Jacob',  expinfo=A002_JK_Condensin_with_plectonemes_expinfo(expi,AllExp(roi));
    case 'Eugene', expinfo=A002_EK_Condensin_with_plectonemes_expinfo(expi,AllExp(roi));
    end    
    
    LoadName=char(strcat(inpath, 'EKMcp_A020_',Exp));             
    disp(strcat('A050_Contentplots: Exps to work through:',num2str(LE-roi)));
    load(strcat(LoadName, '_allresults.mat'));     
    [kymo_duration,kymo_width]=size(kymo_DNA); 
    
    SaveName=[outpath,'EKMcp_A050_roi', num2str(roi)];
    A050_results_roi=content_of_plectonemes(info_DNA,kymo_DNA,SaveName,expinfo);
    A050_results.content_loops_first50(roi)=A050_results_roi.content_loops_first50;
    A050_results.content_loops_last50(roi)=A050_results_roi.content_loops_last50;
    A050_results.tetherlength(roi)=A050_results_roi.tetherlength;
end 
save(strcat(outpath, 'A050_collectedresults.mat'), 'A050_results');
end

load(strcat(outpath, 'A050_collectedresults.mat'), 'A050_results');
if 1
close all;
subplot(1,2,1);
plot(A050_results.tetherlength,A050_results.content_loops_first50, 'o'); hold on;
plot(A050_results.tetherlength,A050_results.content_loops_last50, 'o');
legend('initial', 'final')
 ylabel('total loop %');
 xlabel('tetherlength, pixels');
 ylim([0 100]);
subplot(1,2,2);
plot(A050_results.tetherlength,A050_results.content_loops_last50-A050_results.content_loops_first50, 'o')
 ylabel('final minus first loop %');
 xlabel('tetherlength, pixels');
 ylim([-50 50]);end



  



        