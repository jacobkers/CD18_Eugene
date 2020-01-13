function matrix_arc_sampler
 %JWJK_C:-------------------------------------------------------------------
%Title: %build an arc grid
%Summary: %build arc grid to analyze a side-flow tether
%Input: movie, clicked positions
%Output: re-sampled map with base length & 'arc length' axes
%References: Jacob Kers 2019
%:JWJK_C-------------------------------------------------------------------
close all;
expi=1;  %sim
expi=2; %one-for-all
%expi=3; %181206_150642 RealTimeLoopExtrusion-rot
expi=3; %181125_175911-1
figure(1);
[xa,ya,clickim,source,pickrange]=matrix_getstartpoints(expi,1);
xa
ya
%smooth and equalize;
[xa,ya]=xy_get_smooth_xyline(xa,ya,100,5);
[xa,ya]=xy_equalize_along_contour(xa,ya,200); %best contour

slp=((ya(end)-ya(1))/(xa(end)-xa(1)));       %slope
ybs=(xa-xa(1))*slp+ya(1);                    %baseline

%get a grid of consecutive arcs, each differing just one unit length from
%its former
[unitarcgrid_y,unitarcgrid_x, bestL]=get_unit_length_arcs(xa,ya,ybs,2000);

close all;
subplot(2,2,1);
    pcolor(clickim); colormap hot; shading flat; hold on;
    plot(xa,unitarcgrid_y','w-','MarkerFaceColor', 'w');
    plot(xa,ya,'ro','MarkerFaceColor', 'w');
    title('base curve &unit length sampling curves');
subplot(2,2,2); plot(bestL,'o-');
    ylabel('length');
    xlabel('curve index');
subplot(2,2,4); plot(bestL-round(bestL),'o-');
    ylabel('residu length');
    xlabel('curve index');
  
    
 %% sampling the movie   
info = imfinfo(source);
[ff,~]=size(info);
%sample frames with this grid
[rra,cca]=size(unitarcgrid_x);
hatcurve=prf_build_stetson_mask(1:cca,ceil(cca/6));
hatmask=repmat(hatcurve,1,rra);
for jj=1:ff
    im=flipud(double(imread(source,'Index',jj)));
    filmap=interp2(im,unitarcgrid_x,unitarcgrid_y,'bilinear',0); 
    filmap=matrix_blur_jk(filmap,3);
    filmap=filmap'.*hatmask;   
    prf=sum(filmap(pickrange,:)); prf=prf-min(prf); prf=1*prf/max(prf);
    switch 1
        case 1
        est.x0=40;
        est.b0=0;
        est.N0=1;
        est.psf=8;
        axz=1:length(prf);
        [fit1x1Dfix,YFit]=MLE_One1D_Gaussian_FixPSF(axz,prf,est,0);    
         x0=fit1x1Dfix.x0;
        case 2
        x0=prf_get_periodic_com(prf);
        case 3
            [val,ix]=max(prf);
        x0=ix+prf_subpix_aroundzero(prf);
    end
    pos(jj)=x0;
    figure(2); 
    subplot(2,2,1); pcolor(filmap); colormap hot; shading flat;
    pause(0.05);
    subplot(2,2,3); 
    plot(prf); hold on; 
    %plot(YFit,'r-'); 
    %plot((pos(jj)),prf(round(pos(jj))),'ro'); 
    hold off;
    
    subplot(1,2,2); 
    plot(pos); 
end
outdata=[[1:length(pos)]', pos'];
dlmwrite([pwd,'\arctrack.txt'],outdata,'delimiter','\t');
dum=1;  
    
function [unitarcgrid_y,unitarcgrid_x,bestL]=get_unit_length_arcs(xa,ya,ybs,Nc)
%JWJK_C:-------------------------------------------------------------------
%Summary: %build large amplitude series and pick near-discrete lengths from these
%References: Project Eugene, Jacob Kers 2019
%:JWJK_C-------------------------------------------------------------------

%1) make brute-force curve series
for ii=1:Nc
    amp(ii)=1.25*ii/Nc;
    yxc=ya-ybs;
    yxc2=[yxc(1) ; cumsum(amp(ii)*diff(yxc))+yxc(1)]+ybs;    
    all_CL(ii)=sum(((xa(2:end)-xa(1:end-1)).^2+...
            (yxc2(2:end)-yxc2(1:end-1)).^2).^0.5);  %approximate contour length
end

minL=floor(all_CL(1));
maxL=ceil(all_CL(end));
unitarcgrid_y=zeros(maxL-minL+1,length(xa));
cnt=0;
for LL=minL:maxL
    cnt=cnt+1;
    [val,ix]=min(abs(all_CL-LL));
    unitarcgrid_y(cnt,:)=[[yxc(1) ; cumsum(amp(ix)*diff(yxc))+yxc(1)]+ybs]';
    bestL(cnt)=all_CL(ix);
end
unitarcgrid_x=repmat(xa',cnt,1);
dum=1


function [xa,ya,clickim,source,pickrange]=matrix_getstartpoints(expi,reuseit);
%JWJK_C:-------------------------------------------------------------------
%Title: click or load points
%References: Project Eugene, Jacob Kers 2019
%:JWJK_C-------------------------------------------------------------------
%open sample image; show it
switch expi
    case 1 
        source='C:\Users\jkerssemakers\Dropbox\CD_Data_out\2018_Eugene\2019_09_16 Extrusion simulator\2019_09_16 Extrusion simulator_5.tif'   ;
        pickrange=80:120;
    case 2 
        source='C:\Users\jkerssemakers\CD_Data_in\2018_Eugene\2019_09_18 side_flow_step\one_for_all-1.tif'; 
        pickrange=[60:100];
    case 3 
        source='C:\Users\jkerssemakers\CD_Data_in\2018_Eugene\2019_09_18 side_flow_step\181206_150642 RealTimeLoopExtrusion-rot.tif'; 
        pickrange=110:120;
    case 4 
        source='C:\Users\jkerssemakers\CD_Data_in\2018_Eugene\2019_09_18 side_flow_step\181125_175911-1.tif';
        pickrange=60:200;
end
clickim=flipud(double(imread(source,'Index',20)));

if ~reuseit%click some 10 points;
    pcolor(clickim); colormap hot; shading flat; hold on;
    but=1; cnt=0;
    while but==1
        cnt=cnt+1;
        [x,y,but]=ginput(1);    
        if but~=3, 
            xa(cnt)=x; ya(cnt)=y; 
            plot(xa,ya,'wo-');
        end
    end
else
  [xa,ya]=load_points(expi);  
    
end


function [xa,ya]=load_points(expi);
switch expi
    case 1    
    xa=[24.391       25.847       27.121       28.213       29.851         32.4       35.312       39.681       ...
        45.142       47.872       50.603       52.059       52.969       54.789       55.881];
    ya=[38.197       43.264       48.331       53.859       58.465       64.684       69.981       72.975      ...
        72.745       67.908       63.302       58.695        54.78       50.404       46.719];
    case 2
        xa=[ 11.6037   13.4009   15.1982   17.5945   20.1106   22.9862   25.9816   28.8571   30.0553   30.8940   31.6129   31.7327   32.0922];
        ya=[23.3367   28.7245   33.6633   38.6020   44.4388   48.9286   50.0510   49.1531   45.1122   40.8469   36.1327   32.0918   29.8469];
    case 3
        xa=[43.7108   45.2200   48.8422   50.3514   53.6717   58.5012   62.1233   63.3306   66.0472   68.4620   70.8767   72.6878   75.4044   77.5173   78.7247   79.9320];
        ya=[71.1064   76.8440   82.2230   86.1676   88.3192   87.9606   84.3746   81.5058   76.4854   70.0306   65.0102   59.6312   53.8936   50.3076   47.4388   44.9286]; 
    case 4  
        xa =[34.7258   36.8410   39.4263   43.1866   48.8272   54.7028   59.1682   65.2788   69.5092   72.5645   74.4447   76.3249   77.5000];
        ya =[72.3251   84.8440   92.4213   96.7041   99.0102   96.3746   90.1152   77.5962   65.0773   57.8294   47.6166   41.3571   35.7566];
        
end
