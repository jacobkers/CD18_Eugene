
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
posses_Cnd=PeakFit_kymo(kymo_Cnd,'flatbottom',4);
%posses_DNA=PeakFit_kymo(kymo_DNA,'edgedetect',3);
posses_DNA=PeakFit_kymo(kymo_DNA,'peeling',3);

%% show result
figure(1);
if sho
    [rrc,ccc]=size(kymo_Cnd);
    [rrd,ccd]=size(kymo_DNA);
    subplot(2,2,1); pcolor(kymo_Cnd); shading flat, colormap hot;
    title('Condensin'); ylabel('frame no.');
    subplot(2,2,2); pcolor(kymo_DNA); shading flat, colormap hot;
    title('DNA'); ylabel('frame no.');
    subplot(2,2,3); plot(posses_Cnd(:,2),posses_Cnd(:,1),'ko', 'MarkerSize',2);
    xlabel('position'); ylabel('frame no.');
    xlim([1 ccc]); ylim([1 rrc]);
    subplot(2,2,4); plot(posses_DNA(:,2),posses_DNA(:,1),'ko','MarkerSize',2);
    xlabel('position'); ylabel('frame no.');
    xlim([1 ccd]); ylim([1 rrd]);
end

figure(2);
shiftX=-10;

plot(posses_DNA(:,3)+shiftX,posses_DNA(:,1),'co','MarkerSize',2); hold on;
plot(posses_DNA(:,2)+shiftX,posses_DNA(:,1),'bo','MarkerSize',3); hold on;

%plot(posses_Cnd(:,3),posses_Cnd(:,1),'mo', 'MarkerSize',2); hold on;
plot(posses_Cnd(:,2),posses_Cnd(:,1),'ro', 'MarkerSize',3); hold on;
legend('DNA-1st','DNA-2nd','Condensin-2nd');
xlim([1 ccd]); ylim([1 rrd]);
%% save data
SaveName=char(strcat(generaldatapth, Exp,Expi));
save(strcat(SaveName, '_allresults.mat'),... 
            'kymo_Cnd',     'kymo_DNA',...
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
    
    
function posses=PeakFit_kymo(kymo,peakfitoption, sigma);
posses=[];
[FramesNo,~]=size(kymo);
cleankymo=kymo;
    switch peakfitoption
        case 'peeling'
            for jj=1:FramesNo 
                prf=kymo(jj,:);
                prf=smooth(prf,4);
                [peakprops,buildprf]=peel_peaks_from_profile(prf',2.7,0);   
                cleankymo(jj,:)=buildprf-median(buildprf);
                [betterpeaks, betterpeaksvals]= Refine_Peaks(prf,peakprops(:,3), 0);
                [LL,~]=size(peakprops);    
                posses=[posses; [jj+zeros(LL,1) betterpeaks peakprops(:,3)]];
            end
        case 'flatbottom'
             for jj=1:FramesNo 
                prf=kymo(jj,:);
                 prf=JKD1_PRF_smooth(prf',4);
                [firstpeaks,betterpeaks, betterpeaksvals]=...
                 JKD1_PRF_get1Dpeaksflatbottom(prf,sigma,1,0);%data,sigs,refine,plotit
                dum=1;
                LL=length(betterpeaks);
                posses=[posses; [jj+zeros(LL,1) betterpeaks firstpeaks]];
             end
        case 'edgedetect'
            [kymo_shr, realshrinkfactor]=shrink_kymo(kymo, 1);
          
            [stepfitim_c,stepfitim_r,stepcoords_r,stepcoords_c]=StepfindImage(kymo_shr);
            posses=[stepcoords_c(:,2) realshrinkfactor*stepcoords_c(:,3) realshrinkfactor*stepcoords_c(:,3)];
    end
    
     function [kymo_shr, realshrinkfactor]=shrink_kymo(kymo, shrinkfactor);
     %shorten along x with about the psf in pixels
     
     [rr,cc0]=size(kymo);
     cc1=round(cc0/shrinkfactor);
     kymo_shr=zeros(rr,cc1);
     ax0=1:cc0;
     ax1=linspace(1,cc0,cc1);
     for ii=1:rr
         prf0=kymo(ii,:);
         prf1=interp1(ax0,prf0,ax1,'spline');
         kymo_shr(ii,:)=prf1;
     end
     realshrinkfactor=cc0/cc1;   
         
