
function A020_KymoStackTrack
%JWJK_A:-------------------------------------------------------------------
%Title:A020_KymoStackTrack
%
%Summary: %This function builds a kymograph from images 
%and analyzes it
%JacobKers2018
%
%:JWJK_A-------------------------------------------------------------------

close all;
loadImageJ_kymograph=1;

generaldatapth='D:\jkerssemakers\_Data\CD\2018_Eugene\';  

%% standardized subdirectory names 
Channel_list=[{'DNA\'}, {'Condensin\'}];     %The two subchannels
Kymo_list='kymo_ImageJ\'; %ImageJ-made kymographs
Condensin_Kymo='Kymograph_Condensin.tif';           %if you use it
Dna_Kymo='Kymograph_DNA.tif';                %if you use it';

%% some experiment-specific settings
AllExp=[{'ROI'} {'ROI3'}];  %paths to various experiments   


%% main loop
LE=length(AllExp);  %for all experiments
for ee=1:LE
if mod(ee,2)==0, disp(strcat('Exps to work through:',num2str(LE-ee)));end    
Exp=AllExp(ee);

% 1) build stack of images and kymograph (or load it)
if loadImageJ_kymograph             %load kymograph
    dna_name=char(strcat(generaldatapth, Exp, '\', Kymo_list,Dna_Kymo));
    condensin_name=char(strcat(generaldatapth, Exp,'\', Kymo_list,Condensin_Kymo));
    kymo_DNA=double(imread(dna_name));
    kymo_Cnd=double(imread(condensin_name));    
else                                %make two kymographs
    dna_pth=char(strcat(generaldatapth, Exp,'\', Channel_list(1)));
    condensin_pth=char(strcat(generaldatapth, Exp,'\', Channel_list(2)));    
    kymo_DNA=Build_kymo(dna_pth);
    kymo_Cnd=Build_kymo(condensin_pth);
end

% 2 do peak analysis on each profile of the kymographs
posses_Cnd=PeakFit_kymo(kymo_Cnd,'flatbottom',4);

%2a do loop analysis
kymo_DNA=kymo_DNA-min(kymo_DNA(:));
[kymo_DNA_loop,kymo_DNA_residu]=Clean_Kymo(kymo_DNA);
loopinfo=Analyze_Loop(kymo_DNA_loop,kymo_DNA_residu);

% 3 show result
figure(1);
    [rrc,ccc]=size(kymo_Cnd);
    [rrd,ccd]=size(kymo_DNA);
    subplot(2,2,1); pcolor(kymo_Cnd); shading flat, colormap hot;
    title('Condensin'); ylabel('frame no.');
    subplot(2,2,2); pcolor(kymo_DNA); shading flat, colormap hot;
    title('DNA'); ylabel('frame no.');
    subplot(2,2,3); plot(posses_Cnd(:,2),posses_Cnd(:,1),'ko', 'MarkerSize',2);
    xlabel('position'); ylabel('frame no.');
    xlim([1 ccc]); ylim([1 rrc]);
    %subplot(2,2,4); plot(posses_DNA(:,2),posses_DNA(:,1),'ko','MarkerSize',2);
    %xlabel('position'); ylabel('frame no.');
     %   xlim([1 ccd]); ylim([1 rrd]);

figure(2);
    subplot(1,2,1);
    plot(posses_Cnd(:,2),posses_Cnd(:,1),'ro', 'MarkerSize',6); hold on;
    
    pairx=[loopinfo.left' loopinfo.right'];
    pairy=[loopinfo.fr' loopinfo.fr'];   
    plot(pairx',pairy','b-');
    legend('Condensin','Loop');
    xlim([1 ccd]); ylim([1 rrd]);
    plot(pairx,pairy,'bo','MarkerSize',3,'MarkerFaceColor','b');
    plot(posses_Cnd(:,2),posses_Cnd(:,1),'ro', 'MarkerSize',6); hold on;
    xlabel('position, pixels');
    ylabel('frame no.');
    
    subplot(2,2,2);
    plot(loopinfo.fr, loopinfo.length, 'b-');
    xlabel('frame no.');
    ylabel('looplength, pixels');
     xlim([1 loopinfo.fr(end)]);
    subplot(2,2,4);
    plot(loopinfo.fr, loopinfo.perc, 'k-');
    xlabel('frame no.');
    ylabel('loopcontent, %');
    xlim([1 loopinfo.fr(end)]);

% save data
SaveName=char(strcat(generaldatapth, Exp));
save(strcat(SaveName, '_allresults.mat'),... 
            'kymo_Cnd',     'kymo_DNA',...
            'kymo_DNA_loop', 'kymo_DNA_residu',...
            'posses_Cnd', 'loopinfo');
saveas(gcf,strcat(SaveName, '_plots.jpg'),'jpg');      
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
                %[peakprops,buildprf]=peel_peaks_from_profile(prf',2.7,1);
                [peakprops,buildprf,clusterprops, allclusters]=peel_peaks_from_profile_plusclusters(prf',2.7,1,'NonPeriodic');
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
            [stepfitim_c,stepfitim_r,stepcoords_r,stepcoords_c]=StepfindImage(kymo);
            stepcoords_c=Merge_Edges(stepcoords_c,stepfitim_c, kymo);
            
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
         
function [loop_DNA,residu_DNA]=Clean_Kymo(kymo_DNA);
    %Remove lowerpart(non-condensed) part
    tetherlevel=median(kymo_DNA(:));
    loop_DNA=kymo_DNA-tetherlevel;
    loop_DNA(loop_DNA<0)=0;   %shave off level
    residu_DNA=kymo_DNA-loop_DNA;
    if 0
        close all;
        subplot(2,2,1);
        pcolor(kymo_DNA); shading flat, colormap hot;
        title('original');
        subplot(2,2,2);
        pcolor(loop_DNA); shading flat, colormap hot;
        title('condensed');
        subplot(2,2,3);
        pcolor(residu_DNA); shading flat, colormap hot;
        title('tether residu');
        subplot(2,2,4);
        plot(residu_DNA','k-');  hold on;
        plot(loop_DNA','r-'); 
        title('residu and condensed');
        dum=1;
    end
    
   function loopinfo=Analyze_Loop(loop_kymo,residu_kymo);
       %assume one dna loop ('cluster')
       % 1. set a left and right cluster border, based on treshold from
       % buildcurve
       %2. look for all spots (irrespective of cluster) that fall in this
       %slot. find leftmost and rightmost spot. The are the loop start and
       %stop positons, psf-corrected
       %get out:
       %loopinfo.leftpos: the s
       %loopinfo.rightpos
       %loopinfo.relcontent: amount of signal in the loop compared to the
       %residu
        [FramesNo,~]=size(loop_kymo);
        loopinfo.length=[];
        loopinfo.perc=[];
        loopinfo.fr=[];
        loopinfo.left=[];
        loopinfo.right=[];
        for jj=1:FramesNo 
            loopinfo.fr(jj)=jj;
            prf=loop_kymo(jj,:);
            prf_res=residu_kymo(jj,:);
            prf=(smooth(prf,4));
            [peakprops,buildprf,clusterprops, allclusters]=peel_peaks_from_profile_plusclusters(prf',2.7,0,'NonPeriodic');
            %peakprops: [peakcount PeakVal Xpos Psf ThisSpotFraction(peakcount) CoveredFraction(peakcount) RelChange]];
            %clusterprops: [ii ClusterComPos ClusterContent];
            prf=prf';
            % 1. set a left and right cluster border
            tresh=0.2;
            sel=find(prf>max(prf)*tresh);
            if ~isempty(sel)
                lft_edge=sel(1);
                rgt_edge=sel(end);
                %2 find encompassed spot positions
                inpeaks_ix=find((peakprops(:,3)>lft_edge)&(peakprops(:,3)<rgt_edge));
                inpeaks_x=sort(peakprops(inpeaks_ix,3));
                loopinfo.left(jj)=inpeaks_x(1);
                loopinfo.right(jj)=inpeaks_x(end);
                loopinfo.length(jj)=inpeaks_x(end)-inpeaks_x(1);
                inpeaks_content=nansum(peakprops(inpeaks_ix,5));
                loopinfo.perc(jj)=100*inpeaks_content*sum(prf)/(sum(prf+prf_res));
                end
            dum=1;
        end
     
        dum=1;