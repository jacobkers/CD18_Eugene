function A020_BuildTraces(initval)
% 'Use this section for a Quicksheet'
    %------------------------------------------------------------------
    % This program:   
    %loads image by image
    %loads a n array of spot areas
    %and gets the intensity of these
   %------------------------------------------------------------[JK15]
 % 'End of Quicksheet section'

close all
if nargin<1
    initval=A000_ConfigExp; %Main and image paths per user:
end
spotposses=dlmread(strcat(initval.SaveDataPath,initval.SpotDataResultsName));
N_all=length(spotposses);
imagenames=dir(strcat(initval.ImDataPath,'*.tif*'));
dum=1;
RoiHSize=12;
skips=1;
max_im=10000;
Nims=min(max_im,length(imagenames));;
ff=1:skips:Nims;
nframes=length(ff);
traces_I=zeros(nframes,N_all);
traces_X=zeros(nframes,N_all);
traces_Y=zeros(nframes,N_all);
shiftax=2000*repmat((1:N_all),nframes,1);
for fri=1:nframes
    display(strcat(num2str(nframes-fri+1),' images to go'));
    imidx=ff(fri);
    filnam=imagenames(imidx).name;
    nm=strcat(initval.ImDataPath,filnam);
    FL=double(imread(nm));  %load label- image 
    [rr,cc]=size(FL);
    for spi=1:N_all
        x0=spotposses(spi,1);
        y0=spotposses(spi,2);
        x1=max([1+ RoiHSize, x0]);
        x=min([cc-RoiHSize, x1]);
        y1=max([1+ RoiHSize, y0]);
        y=min([rr-RoiHSize, y1]);
        lox=round(x-RoiHSize);
        hix=round(x+RoiHSize);
        loy=round(y-RoiHSize);
        hiy=round(y+RoiHSize);
        FL_roi=FL(loy:hiy,lox:hix);
%         pcolor(FL_roi); shading flat; pause(0.1);
        spot_fit=SpotsBy_localMaskedROI_iterative(FL_roi);
        traces_I(fri,spi)=spot_fit.N0;
        traces_X(fri,spi)=initval.nmperpix*spot_fit.x0;
        traces_Y(fri,spi)=initval.nmperpix*spot_fit.y0;
    end    
end

%some post_processing
Q1_slot=1:round(nframes/4);
Q4_slot=(round(nframes/4):nframes);
norm_I=mean(traces_I(:));
meanI_norm=mean(traces_I)/norm_I;
ratioI=mean(traces_I(Q4_slot,:))./mean(traces_I(Q1_slot,:));



noflutterspots=traces_I(1,:)>0.1*median(traces_I(1,:));

idx_spots=find(meanI_norm<2 & ratioI<initval.bleachratio_treshold & ratioI>0& noflutterspots);  %'bleachers'
idx_markers=find(meanI_norm>0.5 & ratioI>initval.bleachratio_treshold);

spot_std=(std(traces_X(:,idx_spots)));

if isempty(idx_markers)
    idx_markers=find(meanI_norm>0.3 & ratioI>0.2);
end



traces_X=traces_X-repmat(traces_X(1,:),nframes,1);
traces_Y=traces_Y-repmat(traces_Y(1,:),nframes,1);

%% histograms:
binax=-300:30:300;

XX_spots=traces_X(:,idx_spots); XX_spots=XX_spots(:);
YY_spots=traces_Y(:,idx_spots); YY_spots=YY_spots(:);
XX_markers=traces_X(:,idx_markers); XX_markers=XX_markers(:);
YY_markers=traces_Y(:,idx_markers); YY_markers=YY_markers(:);

N_spots=length(idx_spots);
N_markers=length(idx_markers);
N_meas_spots=length(XX_spots);
N_meas_markers=length(XX_markers);
histX_spots=hist(XX_spots, binax)/N_meas_spots*100;
histY_spots=hist(YY_spots, binax)/N_meas_spots*100;
histX_markers=hist(XX_markers, binax)/N_meas_markers*100;
histY_markers=hist(YY_markers, binax)/N_meas_markers*100;

[FWHM_spots_X]=get_FWHM_mainpeak(histX_spots,binax);
[FWHM_markers_X]=get_FWHM_mainpeak(histX_markers,binax);

%% plotting
figure(101);
subplot(1,2,1); 
    plot(traces_X(:,idx_spots), 'r-');    hold on; 
    plot(traces_X(:,idx_markers), 'b-');    hold on; 
    title(Replace_underscores([initval.BleachCurveResultsName 'traces']));
    xlabel('frame');
    ylabel('X-pos, nm');
subplot(1,2,2);
    plot(meanI_norm,ratioI, 'y*'); hold on;
    plot(meanI_norm(idx_spots),ratioI(idx_spots), 'ro');
    plot(meanI_norm(idx_markers),ratioI(idx_markers), 'bo');
    title('trends');
    ylim([-1 2.5]);
    legend('all', 'spots', 'markers')
    xlabel('mean normalized intensity, a.u');
    ylabel('bleach_ratio, a.u.');

saveas(gcf, [initval.SaveDataPath, 'A020_', initval.BleachCurveResultsName 'trends.jpg']);    

figure(93) ; 
subplot(2,2,1); 
    plot(traces_I(:,idx_spots), 'r-');    hold on; 
    plot(traces_I(:,idx_markers), 'b-');    hold on; 
    title(Replace_underscores([initval.BleachCurveResultsName 'traces']));
    text(0.6*nframes,max(max(traces_I(:,idx_markers))),[num2str(N_spots) ' spots']);
    text(0.6*nframes,0.8*max(max(traces_I(:,idx_markers))),[num2str(N_markers) ' markers']);
    xlabel('frame');
    ylabel('fluorescent intensity, a.u.');
subplot(2,2,2); 
    plot(traces_X(:,idx_spots), traces_Y(:,idx_spots), 'ro');    hold on;
    plot(traces_X(:,idx_markers), traces_Y(:,idx_markers),'bo');    hold on;
    pause(0.05);
    title('positions');
    xlabel('x-position, nm');
    ylabel('y-position, nm');
    axis equal; 
 subplot(2,2,3); 
    bar(binax, histX_spots, 'r'); hold on; 
    plot(binax, histX_markers, 'b-');
    title('X-position');
    ylabel('frequency, %');
    xlabel('position, nm');
 subplot(2,2,4);
    bar(binax, histY_spots, 'r'); hold on; 
    plot(binax, histY_markers, 'b-');
    title('Y-position');
    ylabel('frequency, %');
    xlabel('position, nm');

saveas(gcf, [initval.SaveDataPath, 'A020_', initval.BleachCurveResultsName 'overview.jpg']);

%% for paper
%1) rebuttal panel
%2) make named data for saving to .mat:
traces_spots_X=traces_X(:,idx_spots);
traces_spots_Y=traces_Y(:,idx_spots);
traces_spots_frames=repmat((1:nframes)',1,length(idx_spots)); 
traces_spots_IDs=repmat((1:length(idx_spots)),nframes,1); 
traces_markers_X=traces_X(:,idx_markers); 
traces_markers_Y=traces_Y(:,idx_markers);
traces_markers_frames=repmat((1:nframes)',1,length(idx_markers)); 
traces_markers_IDs=repmat((1:length(idx_markers)),nframes,1); 

figure(136);
subplot(1,2,1); 
    plot(traces_spots_X,    traces_spots_Y, 'ro');    hold on;
    plot(traces_markers_X,  traces_markers_Y,'bo');    hold on;
    pause(0.05);
    title('positions');
    xlabel('x-position, nm');
    ylabel('y-position, nm');
    axis equal; 
 subplot(2,2,2); 
    bar(binax, histX_spots, 'r'); hold on; 
    title('ori motion');
    ylabel('frequency, %');
    xlabel('position, nm');
    text(200,0.5*max(histX_spots),'FWHM:')
    text(200,0.3*max(histX_spots),[num2str(round(FWHM_spots_X)), 'nm']);
    
subplot(2,2,4); 
    bar(binax, histX_markers, 'b'); hold on; 
    title('cell motion');
    ylabel('frequency, %');
    xlabel('position, nm');
    text(200,0.5*max(histX_spots),'FWHM:')
    text(200,0.3*max(histX_spots),[num2str(round(FWHM_markers_X)), 'nm']);
saveas(gcf, [initval.SaveDataPath, 'rebuttal_23\A020_', initval.BleachCurveResultsName '_rebuttal_spotdynamics.jpg']);
save([initval.SaveDataPath, 'rebuttal_23\A020_', initval.BleachCurveResultsName '_rebuttal_spotdynamics.mat'],...
    'traces_spots_X' ,'traces_spots_Y', ...
    'traces_markers_X', 'traces_markers_Y',...
    'binax', 'histX_spots', 'histX_markers')

%%excell data:
%spots:
header1=[{'frame'},{'spot_ID'},{'traces_spots_X'},{'traces_spots_Y'}];
data1=[traces_spots_frames(:), traces_spots_IDs(:), traces_spots_X(:), traces_spots_Y(:)];
xlswrite([initval.SaveDataPath, 'rebuttal_23\A020_', initval.BleachCurveResultsName '_rebuttal_spotdynamics.xls'], header1, 'spots', 'A1');
xlswrite([initval.SaveDataPath, 'rebuttal_23\A020_', initval.BleachCurveResultsName '_rebuttal_spotdynamics.xls'], data1, 'spots', 'A2');
%markers:
header2=[{'frame'},{'spot_ID'},{'traces_markers_X'},{'traces_markers_Y'}];
data2=[traces_markers_frames(:), traces_markers_IDs(:), traces_markers_X(:), traces_markers_Y(:)];
xlswrite([initval.SaveDataPath, 'rebuttal_23\A020_', initval.BleachCurveResultsName '_rebuttal_spotdynamics.xls'], header2, 'markers', 'A1');
xlswrite([initval.SaveDataPath, 'rebuttal_23\A020_', initval.BleachCurveResultsName '_rebuttal_spotdynamics.xls'], data2, 'markers', 'A2');
%histograms
header3=[{'bin'},{'spots, percentage'},{'markers,  percentage'}];
data3=[binax', histX_spots',histX_markers'];
xlswrite([initval.SaveDataPath, 'rebuttal_23\A020_', initval.BleachCurveResultsName '_rebuttal_spotdynamics.xls'], header3, 'histograms', 'A1');
xlswrite([initval.SaveDataPath, 'rebuttal_23\A020_', initval.BleachCurveResultsName '_rebuttal_spotdynamics.xls'], data3, 'histograms', 'A2');



dlmwrite(strcat(initval.SaveDataPath,[initval.BleachCurveResultsName 'bleachcurve.txt']),traces_I);
dlmwrite(strcat(initval.SaveDataPath,[initval.BleachCurveResultsName 'X.txt']),traces_X);
dlmwrite(strcat(initval.SaveDataPath,[initval.BleachCurveResultsName 'Y.txt']),traces_Y);
disp(strcat(num2str(N_all), ' traces made'));



