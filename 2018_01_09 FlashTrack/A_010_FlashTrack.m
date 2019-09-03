function A_010_FlashTrack
%JWJK_A:-------------------------------------------------------------------
%Title: FlashTrack
%Project: CD lab, researchers Mahipal, Eugene. Written by: Jacob
%Summary: %This function does quick tracking of images containing one spot
%Approach: Cut Rois by pre-estimates, clip circular mask to measure spot
%content
%Input and how to run it: fill directories with movie frame images of 
%individual spots;
%Output: table with x and y coordinates, graphs. 
%References: MH, Science 360, 102-105 (2018)
%:JWJK_A-------------------------------------------------------------------


close all;
sho=1;  %intermdiate showing (slows sown)
skip=1; %just for quick runs otherwise leave at 1
wd=2;  %median filter width
CropSeries=1E6;

codepth=pwd;
savepth=strcat(codepth,'\GeneralResults\');

Expi='EK_Develop multispot';

switch Expi
    case 'Test', %test
        datapth=codepth;
        SourceList={strcat(datapth,'ROI_for_Tracking\')};
    case 'Trace180107_181053_C2'  %test
        datapth='K:\bn\cd\Shared\SMF_MG_JK_AK\Condensin\180107 _HILO_Atto647N_Condensin\180107_181053_Analysis\DNA_Condensin_2\';
        %datapth='K:\bn\cd\Shared\SMF_MG_JK_AK\Condensin\180107 _HILO_Atto647N_Condensin\180107_184735_Analysis\DNA_condensin_4\';
        SourceList={strcat(datapth,'Condensin_ROI\')};
    case 'FullSeries'
        SourceList=SourcePathList;  %data mining on K!
        LE=length(SourceList);
        if CropSeries<LE, SourceList=SourceList(1:CropSeries);end
    case 'EK_Develop multispot'
         datapth='D:\jkerssemakers\_Recent\CD\BN_CD18_Eugene\2018_01_09 FlashTrack\';        
        SourceList={strcat(datapth,'Condensin_ROI\')};
end

LE=length(SourceList);
for ee=1:LE
if mod(ee,2)==0, disp(strcat('Exps to work through:',num2str(LE-ee),'with every',num2str(skip),'frame skipped'));end    

sourcepth=char(SourceList(ee));
cd(sourcepth); imlist=dir('*.tif'); cd(codepth);
FramesNo=length(imlist);

%% 1) build stack of images and averaged image
for jj=1:skip:FramesNo
    if mod(jj,50)==0, disp(strcat('Stack frames to load:',num2str(FramesNo-jj)));end
    source=strcat(sourcepth,imlist(jj).name);
    pic=(double(imread(source)));
    if jj==1
         [rr,cc]=size(double(pic));
         ImStack=zeros(rr,cc,FramesNo);               
    end
    ImStack(:,:,jj)=pic;
    
end

%3 track center of roi 
xtrace=zeros(ceil(FramesNo/skip),1);
ytrace=zeros(ceil(FramesNo/skip),1);
Counts=zeros(ceil(FramesNo/skip),1);
for jj=1:skip:FramesNo      %for all frames
    if mod(jj,50)==0, disp(strcat('Track frames to go:',num2str(FramesNo-jj)));end
    roi0=ImStack(:,:,jj);  
   
    %pre-find max
    prfx=smooth(sum(roi0),3); 
    prfy=smooth(sum(roi0'),3); 
    
    %Build edge blocker mask
    if jj==1
        [rr,cc]=size(roi0);
        XStetsonCurve=Build1DStetsonMask(prfx);
        YStetsonCurve=Build1DStetsonMask(prfy);
        XYStetsonMap=repmat(XStetsonCurve',rr,1).*repmat(YStetsonCurve,1,cc);
    end
    
    %block edges
    prfx=(prfx-median(prfx)).*XStetsonCurve;
    prfy=(prfy-median(prfy)).*YStetsonCurve;
    roi0=(roi0-median(roi0(:))).*XYStetsonMap;
    
    [~,yc]=max(prfy);
    [~,xc]=max(prfx);
    
    [roi,xoff,yoff]=Cut_roi(roi0,xc,yc);
    [x1,y1,~,~]=TrackXY_by_2DXCor(roi);

    [~,~,Ispot,Ibackground_level,spotim_clipped,spotim_bc,bckim]=DoubleMaskedCom(roi,x1,y1);

    idx=ceil(jj/skip);
    xtrace(idx,:)=x1+xoff;
    ytrace(idx,:)=y1+yoff; 
    
    Counts(idx)=Ispot;
    if sho
    subplot(2,3,1); pcolor(roi0); colormap hot; shading flat; axis equal; axis tight; hold on;
    plot(x1+xoff+0.5,y1+yoff+0.5,'bo-','MarkerSize',6);  hold off
    title('original');    
    subplot(2,3,2); pcolor(roi); colormap hot; shading flat; axis equal; axis tight; 
    title('cut');
    subplot(2,3,3); pcolor(spotim_bc); colormap hot; shading flat; axis equal; axis tight; hold on;
    plot(x1+0.5,y1+0.5,'bo-','MarkerSize',6);  hold off
    title('clipped') 
    
    subplot(2,2,3); 
    plot(xtrace(1:jj),'-'); hold on; 
    plot(ytrace(1:jj),'r-'); hold off;
    xlabel('frame no.');
    ylabel('pos, pixels');
    legend('X' , 'Y');
    subplot(2,2,4); plot(Counts(1:jj),'k-');
    xlabel('frame no.');
    ylabel('counts');
    pause(0.1);
    end
end


xtrace=medfilt1(xtrace,wd); 
ytrace =medfilt1(ytrace,wd); 
Counts=medfilt1(Counts,wd);


figure(2);
subplot(2,1,1); 
plot(xtrace,'-'); hold on;
plot(ytrace,'r-'); hold off;
    xlabel('frame no.');
    ylabel('pos, pixels');
    legend('X' , 'Y');
subplot(2,1,2); plot(Counts,'k-');
    xlabel('frame no.');
    ylabel('counts');
    pause(0.1);
    %[~]=ginput(1);

    %cleaning     
xtrace=xtrace-xtrace(1);
ytrace=ytrace-ytrace(1);

%% saving
TraceName=strcat('TrackResults');

dlmwrite(strcat(sourcepth,TraceName,'.txt'),[xtrace ytrace Counts]);
saveas(gcf,strcat(sourcepth,TraceName,'jpg'), 'jpg');

%Counts
PL=ceil(sqrt(LE));

figure(3); subplot(PL,PL,ee); plot(Counts,'k-'); hold on; text(length(Counts),0,num2str(ee),'FontSize',6);
set(gca,'FontName','Arial','FontSize',6,'FontUnits', ...
'points');


%X
figure(4); subplot(PL,PL,ee); plot(xtrace,'b-'); hold on; text(length(Counts),0,num2str(ee),'FontSize',6);
set(gca,'FontName','Arial','FontSize',6,'FontUnits', ...
'points');

%Y
figure(5); subplot(PL,PL,ee); plot(ytrace,'r-'); hold on; text(length(Counts),0,num2str(ee),'FontSize',6);
set(gca,'FontName','Arial','FontSize',6,'FontUnits', ...
'points');

end
saveas(figure(3),strcat(savepth,Expi,'_Alltraces_Counts.jpg'), 'jpg');
saveas(figure(3),strcat(savepth,Expi,'_Alltraces_Counts'));
saveas(figure(4),strcat(savepth,Expi,'_Alltraces_Xposition.jpg'), 'jpg');
saveas(figure(4),strcat(savepth,Expi,'_Alltraces_Xposition'));
saveas(figure(5),strcat(savepth,Expi,'_Alltraces_Yposition.jpg'), 'jpg');
saveas(figure(5),strcat(savepth,Expi,'_Alltraces_Yposition'));

function SourceList=SourcePathList
SourceList=[...
{'K:\bn\cd\Shared\SMF_MG_JK_AK\Matlab\2018_01_09 FlashTrack\ROI_for_Tracking\'};                                                    %1
{'K:\bn\cd\Shared\SMF_MG_JK_AK\Condensin\180107 _HILO_Atto647N_Condensin\180107_184735_Analysis\DNA_condensin_4\'};                 %2
{'K:\bn\cd\Shared\SMF_MG_JK_AK\Condensin\180107 _HILO_Atto647N_Condensin\180107_165026 ananlysis\DNA_condensin_1\Condensin_ROI\'};  %3
{'K:\bn\cd\Shared\SMF_MG_JK_AK\Condensin\180107 _HILO_Atto647N_Condensin\180107_165026 ananlysis\DNA_condensin_2\Condensin_ROI\'};  %4
{'K:\bn\cd\Shared\SMF_MG_JK_AK\Condensin\180107 _HILO_Atto647N_Condensin\180107_165026 ananlysis\DNA_condensin_3\Condensin_ROI\'};  %5
{'K:\bn\cd\Shared\SMF_MG_JK_AK\Condensin\180107 _HILO_Atto647N_Condensin\180107_165026 ananlysis\DNA_condensin_4\Condensin_ROI\'};  %6
{'K:\bn\cd\Shared\SMF_MG_JK_AK\Condensin\180107 _HILO_Atto647N_Condensin\180107_165817_Analysis\DNA_Condensin_1\Condensin_ROI\'};   %7
{'K:\bn\cd\Shared\SMF_MG_JK_AK\Condensin\180107 _HILO_Atto647N_Condensin\180107_165817_Analysis\DNA_Condensin_2\Condensin_ROI\'};   %8
{'K:\bn\cd\Shared\SMF_MG_JK_AK\Condensin\180107 _HILO_Atto647N_Condensin\180107_180147 Analysis\DNA_Condensin_1\Condensin_ROI\'};   %9
{'K:\bn\cd\Shared\SMF_MG_JK_AK\Condensin\180107 _HILO_Atto647N_Condensin\180107_181053_Analysis\DNA_Condensin_1\Condensin_ROI_2\'}; %10
{'K:\bn\cd\Shared\SMF_MG_JK_AK\Condensin\180107 _HILO_Atto647N_Condensin\180107_181053_Analysis\DNA_Condensin_1\Condensin_ROI\'};   %11
{'K:\bn\cd\Shared\SMF_MG_JK_AK\Condensin\180107 _HILO_Atto647N_Condensin\180107_181053_Analysis\DNA_Condensin_2\Condensin_ROI\'};   %12
{'K:\bn\cd\Shared\SMF_MG_JK_AK\Condensin\180107 _HILO_Atto647N_Condensin\180107_181053_Analysis\DNA_Condensin_3\Condensin_ROI\'};   %13
{'K:\bn\cd\Shared\SMF_MG_JK_AK\Condensin\180107 _HILO_Atto647N_Condensin\180107_181053_Analysis\DNA_Condensin_4\Condensin_ROI\'};   %14
{'K:\bn\cd\Shared\SMF_MG_JK_AK\Condensin\180107 _HILO_Atto647N_Condensin\180107_182044_Analysis\DNA_Condensin_1\Condensin_ROI\'};   %15
{'K:\bn\cd\Shared\SMF_MG_JK_AK\Condensin\180107 _HILO_Atto647N_Condensin\180107_182044_Analysis\DNA_Condensin_2\Condensin_ROI\'};   %16
{'K:\bn\cd\Shared\SMF_MG_JK_AK\Condensin\180107 _HILO_Atto647N_Condensin\180107_182044_Analysis\DNA_Condensin_3\Condensin_ROI\'};   %17
{'K:\bn\cd\Shared\SMF_MG_JK_AK\Condensin\180107 _HILO_Atto647N_Condensin\180107_182044_Analysis\DNA_Condensin_4\Condensin_ROI_1\'}; %18
{'K:\bn\cd\Shared\SMF_MG_JK_AK\Condensin\180107 _HILO_Atto647N_Condensin\180107_184735_Analysis\DNA_Condensin_2\Condensin_ROI\'};   %19
{'K:\bn\cd\Shared\SMF_MG_JK_AK\Condensin\180107 _HILO_Atto647N_Condensin\180107_184735_Analysis\DNA_Condensin_4\Condensin_ROI\'};   %20
{'K:\bn\cd\Shared\SMF_MG_JK_AK\Condensin\180107 _HILO_Atto647N_Condensin\180107_184735_Analysis\DNA_Condensin_5\Condensin_ROI\'};   %21
{'K:\bn\cd\Shared\SMF_MG_JK_AK\Condensin\180107 _HILO_Atto647N_Condensin\180107_184735_Analysis\DNA_Condensin_6\Condensin_ROI_1\'}; %22
{'K:\bn\cd\Shared\SMF_MG_JK_AK\Condensin\180107 _HILO_Atto647N_Condensin\180107_190027_Analysis\DNA_condensin_2\Condensin_ROI\'};   %23
{'K:\bn\cd\Shared\SMF_MG_JK_AK\Condensin\180107 _HILO_Atto647N_Condensin\180107_190027_Analysis\DNA_condensin_3\Condensin_ROI\'};   %24
{'K:\bn\cd\Shared\SMF_MG_JK_AK\Condensin\180107 _HILO_Atto647N_Condensin\180107_191013_Analysis\DNA_Condensin_Analysis\Condensin_ROI\'};%25
{'K:\bn\cd\Shared\SMF_MG_JK_AK\Condensin\180107 _HILO_Atto647N_Condensin\180107_190027_Analysis\DNA_Condensin_1\Condensin_ROI1\'};  %26
{'K:\bn\cd\Shared\SMF_MG_JK_AK\Condensin\180107 _HILO_Atto647N_Condensin\180107_165026 ananlysis\DNA_condensin_4\Condensin_ROI\'};  %27
{'K:\bn\cd\Shared\SMF_MG_JK_AK\Condensin\180107 _HILO_Atto647N_Condensin\180107_165026 ananlysis\BareDNA for condensin binding time analysis\Condensin_ROI\correct\'}; %28
{'K:\bn\cd\Shared\SMF_MG_JK_AK\Condensin\180107 _HILO_Atto647N_Condensin\180107\180107_212400_Analysis\DNA_condensin_1\Condensin_ROI\'}; %29
{'K:\bn\cd\Shared\SMF_MG_JK_AK\Condensin\180107 _HILO_Atto647N_Condensin\180107_181053_Analysis\DNA_Condensin_2\Condensin_ROI_unbinding\'}; %30
{'K:\bn\cd\Shared\SMF_MG_JK_AK\Condensin\180107 _HILO_Atto647N_Condensin\180107\180107_233307\DNA_Binding_Condensin\'};             %31
];


