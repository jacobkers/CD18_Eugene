
function A020_KymoStackTrack
%JWJK_A:-------------------------------------------------------------------
%Title:A020_KymoStackTrack
%
%Summary: %This function builds a kymograph from images 
%and analyzes it
%JacobKers2018
%
%Approach: 
%
%Input
%
%Output
%
%References
%
%:JWJK_A-------------------------------------------------------------------

close all;
sho=1;  %intermdiate showing (slows sown)
wd=2;  %median filter width
CropSeries=1E6;
codepth=pwd;
loadImageJ_kymograph=1;

generaldatapth='D:\jkerssemakers\_Recent\CD\BN_CD18_Eugene\';        
AllExp=[{'ROI\'}];  %paths to various experiments      

Expi='Test_EK01';

switch Expi
    case 'Test_EK01', %test
         
         Channel_list=[{'DNAROI\'}, {'CondensinROI\'}];     %The two channels
         Kymo_list='kymo_ImageJ\'; %ImageJ-made kymographs
         Condensin_Kymo='Kymograph_Condensin.tif';           %if you use it
         Dna_Kymo='Kymograph_DNA.tif';                %if you use it';
end

LE=length(AllExp);  %for all experiments
for ee=1:LE
if mod(ee,2)==0, disp(strcat('Exps to work through:',num2str(LE-ee),'with every',num2str(skip),'frame skipped'));end    
Exp=AllExp(ee);

%% 1) build stack of images and kymograph (or load it)
if loadImageJ_kymograph             %load kymograph
    dna_name=char(strcat(generaldatapth, Exp, Kymo_list,Dna_Kymo));
    condensin_name=char(strcat(generaldatapth, Exp, Kymo_list,Condensin_Kymo));
    kymo_DNA=double(imread(dna_name));
    kymo_Cnd=double(imread(condensin_name));    
else                                %make two kymographs
    dna_pth=char(strcat(generaldatapth, Exp, Channel_list(1)));
    condensin_pth=char(strcat(generaldatapth, Exp, Channel_list(2)));    
    kymo_DNA=Build_kymo(dna_pth);
    kymo_Cnd=Build_kymo(condensin_pth);
end

%% 2 do peak analysis on each profile of the kymographs
[posses_DNA,kymo_DNA_cln]=PeakFit_kymo(kymo_DNA);
[posses_Cnd,kymo_Cnd_cln]=PeakFit_kymo(kymo_Cnd);

%% show result
figure(1);
if sho
    subplot(2,2,1); pcolor(kymo_Cnd); shading flat, colormap hot;
    subplot(2,2,2); pcolor(kymo_DNA); shading flat, colormap hot;
    subplot(2,2,3); pcolor(kymo_Cnd_cln); shading flat, colormap hot; hold on;
    subplot(2,2,4); pcolor(kymo_DNA_cln); shading flat, colormap hot; hold on;    
    subplot(2,2,3); plot(posses_Cnd(:,2),posses_Cnd(:,1),'wx', 'MarkerSize',1);
    subplot(2,2,4); plot(posses_DNA(:,2),posses_DNA(:,1),'wx','MarkerSize',1);
end

%% save data
SaveName=char(strcat(generaldatapth, Exp,Expi));
save(strcat(SaveName, '_allresults.mat'),... 
            'kymo_Cnd',     'kymo_DNA',...
            'kymo_Cnd_cln','kymo_DNA_cln',...
            'posses_Cnd','posses_DNA');
end

function kymo=Build_kymo(pth)
    codepth=pwd;
    cd(pth); imlist=dir('*.tif'); cd(codepth);
    FramesNo=length(imlist);
    for jj=1:FramesNo  %for all frames        
        if mod(jj,50)==0, disp(strcat('Stack frames to load:',num2str(FramesNo-jj)));end
        source=strcat(pth,imlist(jj).name);
        pic=(double(imread(source)));
        if jj==1
             [rr,cc]=size(double(pic));
             kymo=zeros(FramesNo,rr);         
        end
        prf=sum(pic');
        prf=prf-min(prf);
        kymo(jj,:)=prf;   
    end
    
    
function [posses,cleankymo]=PeakFit_kymo(kymo);
posses=[];
[FramesNo,~]=size(kymo);
cleankymo=kymo;
for jj=1:FramesNo
    prf=kymo(jj,:);
    [peakprops,buildprf]=peel_peaks_from_profile(prf,2.4,0);   
    cleankymo(jj,:)=buildprf-median(buildprf);
    [LL,~]=size(peakprops);
    posses=[posses; [jj+zeros(LL,1) peakprops(:,3)]];
end
