function Psf_meas=B050_EstimatePSF(batchrunindex)
%JWJK_B:-------------------------------------------------------------------
%Estimate Effective pointspread functionfor  cluster analysis
%
%Summary: This function estimates the effective point spread function for 
%a set of images.
%
%Approach: 
%
%Input: data in .mat files stored in former analysis steps.
%
%Output: estimate of psf.
%
%
%:JWJK_B-------------------------------------------------------------------


close all;
savit=1;


if nargin<1 
    batchrunindex=1;
end

initval=A000__WF_Get_JacobPathsandExperiments(batchrunindex);
allframes=length(initval.Cell_Labels);

goodcount=0;
ClusterMaxNo=0;
MatFilePath=strcat(initval.resultpath,'ResultsPerCellMatlab',initval.DirSep);
all_d_peakslopes=[]; all_d_slopeslopes=[];
disp(strcat('PSF_Analysis..', num2str(allframes), 'cells to go'));
for jj=1:allframes    
    cellno=char(initval.Cell_Labels{jj});     
    CellName=strcat('ResultsOfCell',cellno,'.mat'); 
        load(strcat(MatFilePath,CellName));   
        if GeneralCellProps.Okayproduct|initval.PassAllCells
            %disp('good cell')
%            disp(strcat('PSF_Analysis..', num2str(allframes-jj+1), 'cells to go'));
            goodcount=goodcount+1;        
            cellno=char(initval.Cell_Labels{jj});    
            NumCellLabel=BuildNumericCellLabel(cellno);
            CellSpecs=[goodcount NumCellLabel];      
            chro_pic=Remove_Background(chro_pic,'Min');  
            Psf_est=initval.Psf_est;
            
            if 0 %strcmp(cellno,'100233'),
                sho=1;
                close(gcf);
                pcolor(chro_pic); shading flat; colormap bone; hold on;
                axis equal
                [~]=ginput(1);
            else sho=0;
            end
            
            d_peakslopes=Measure_SlopeDistances(chro_pic,Psf_est,sho);
            all_d_peakslopes=[all_d_peakslopes; d_peakslopes];

        end        
end

%all_d_peakslopes=all_d_peakslopes*2/pi;
%Apply correction since most edges will be crossed at oblique angles;
%correction is just cosine integration(tends to underestimate psf)

LPS=length(all_d_peakslopes);
scatax=1.0*rand(LPS,1)-0.4;

slopedata=all_d_peakslopes+scatax;


Psf_peak=nanmedian(slopedata);


 hx=0:0.5:10;
 Histdist=hist(slopedata,hx);
[flag,cleandata]=JKD1_PRF_outlier_flag(slopedata,3,0.7,'positive',0);
OneSigma=std(cleandata);
Psf_meas=Psf_peak-OneSigma;  %steepest cutoff



if  1% nargin <1
    subplot(2,1,1);
    plot(all_d_peakslopes+scatax,'bo'); pause(0.02); hold on;
    xlabel('index,a.u.)')
    ylabel('distance+scatter, pixel units');
    ylim([0 12]);
    
    subplot(2,1,2);
    bar(hx,Histdist,'b');
    xlabel('distance, pixel units)')
    ylabel('counts');
    text(6,0.8*max(Histdist),['PSF:',num2str(Psf_meas)]);
    pause(0.5); %[~]=ginput(1);
end



function d_peakslopes=Measure_SlopeDistances(pic,Psf0,sho)
%JWJK_C:-------------------------------------------------------------------
%Summary: Get distances between peaks and valleys and steepest slopes
%
%Approach: 
%
%Input: data in .mat files stored in former analysis steps.
%
%Output: estimate of psf.
%
%
%:JWJK_C-------------------------------------------------------------------
%smooth noise out via estimate
    blur_pic=JKD2_IM_smoothJK(pic,ceil(Psf0/3));
    blur_pic=blur_pic-min(blur_pic(:));
    subplot(1,2,1);
    d_peakslopesx=Get_PeakMaxSlopeDistances(blur_pic,sho);
    subplot(1,2,2);
    d_peakslopesy=Get_PeakMaxSlopeDistances(blur_pic',sho);   
    d_peakslopes=[d_peakslopesx; d_peakslopesy]; 
    
    dum=1;
    
    
    
    
function d_peakslopes=Get_PeakMaxSlopeDistances(blur_pic,sho);
    %find inflection distances along rows of picture: first find peaks,
    %then look for nearby up- and downward slopes
        
    data=blur_pic(:);
    data_dif=diff(data);
    pks=Get_Maxes(data);  %indices of local maxes 1D    
    maxslopes=Get_Maxes(data_dif)+1;  %indices of local maxes 1D
    minslopes=Get_Maxes(-data_dif)+1;  %indices of local minima 1D 
    [dupletsUp,dupletsDown]=Get_nearest_inflectiondistances(pks,maxslopes,minslopes,data);    
    d_peakslopes=[(dupletsUp(:,2)-dupletsUp(:,1));...
                  (dupletsDown(:,2)-dupletsDown(:,1))];
   if sho
       [rr,cc]=size(blur_pic);
       [XX,YY]=meshgrid(1:cc,1:rr);  %coordinates
       
       pcolor(blur_pic'); shading flat; colormap bone; hold on;
       upbarsX=XX(dupletsUp);
       upbarsY=YY(dupletsUp);
       plot(upbarsY',upbarsX','y-','Linewidth',1);
       downbarsX=XX(dupletsDown);
       downbarsY=YY(dupletsDown);
       plot(downbarsY',downbarsX','y-','Linewidth',1);
       
       plot(downbarsY(:,1),downbarsX(:,1),'yo','MarkerFacecolor','k');
       plot(upbarsY(:,2),upbarsX(:,2),'yo','MarkerFacecolor','k')
       axis equal; axis off;
       dum=1;
   end
    
    
 function [dupletsUp,dupletsDown]=Get_nearest_inflectiondistances(pks,maxslopes,minslopes,data);
    %output: upwardslope, pk, downwardslope (if accepted)
    LE=length(pks);  %use only peaks
    valmax=max(data);
    dupletsUp=zeros(LE,2);
    dupletsDown=zeros(LE,2);
    for ii=1:LE   
        leftupslopes=maxslopes(maxslopes<pks(ii));
        [distMx,idxMx]=min(pks(ii)-leftupslopes);  %nearest upward slope on left
        rightdownslopes=minslopes(minslopes>pks(ii));
        [distMn,idxMn]=min(-pks(ii));  %nearest downward slope on right
        dupletsUp(ii,:)=[leftupslopes(idxMx) pks(ii)];       %up-peak   
        dupletsDown(ii,:)=[pks(ii) rightdownslopes(idxMn)];  %peak-down  
    end    
        %clean up from too shallow peaks
        %upward slopes
        minrange=0.25;  %of max
        
        valslopes=data(dupletsUp(:,1));  %value of local peak
        valpeaks=data(dupletsUp(:,2));     %value of local slope
        valdifs=abs(valslopes-valpeaks);       %difference        
        dupletsUp=dupletsUp((valdifs>minrange*valmax),:);  %remove too shallow ones
        
        %downward slopes
        valslopes=data(dupletsDown(:,2));  %value of local peak
        valpeaks=data(dupletsDown(:,1));     %value of local slope
        valdifs=abs(valslopes-valpeaks);       %difference        
        dupletsDown=dupletsDown((valdifs>minrange*valmax),:);  %remove too shallow ones
        dum=1;
        

function idxes=Get_Maxes(data);
            idxes=find(data(2:end-1)>data(1:end-2) &...
             data(2:end-1)>data(3:end))+1;  %indices of local maxes 1D
%            plot(data(idxes),'o');
%             dum=1;
