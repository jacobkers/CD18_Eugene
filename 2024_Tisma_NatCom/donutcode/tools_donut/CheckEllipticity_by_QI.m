    
function QI=CheckEllipticity_by_QI(im,QI, sho);
%This function : first performs a QI centering routine;
%Then evaluates the final radial pattern by expansion-fitting each with its own
%radial average. %The angle-dependent expansion factor tells us how large
%elliptic deviations are.
close all
 EllipticityVsFocus=[];
 cnt=0;
cnt=cnt+1;
  if nargin <3
      sho=1
      pattern='RealEggy'; 'Sim' %;'Sim'; 'RealEggy' ; %'Sim'; %'RealEggy';
      switch pattern
          case 'Sim'
            test.PicSize=100; 
            x0=test.PicSize/2; 
            y0=test.PicSize/2; 
            test.PatternRingRadius=test.PicSize/8;
            test.PatternRingWidth=test.PicSize/30;
            im=MakeSemiCircularPattern(x0,y0,test,'Flower'); %'Eggshaped'; 'Elliptoid' ; 'Flower'
          case  'RealEggy'
            impth='D:\jkerssemakers\My Documents\BN ND Data All\Tracking_Jordi\RealBeadsJordi\Bead75\';
            files = dir(strcat(impth,'*.jpg')); 	
            [k,l]=size(files);   
            imname=files(800).name; %300
            im=double(imread(strcat(impth,imname))); 
      end
   
    [rr,cc]=size(im);
    QI.radialoversampling=2;
    QI.angularoversampling=0.7;
    QI.minradius=0;
    QI.maxradius=rr/3;
    QI.iterations=5;
    QI=TrackXY_by_QI_Init(im);
  end  
 %close all   

 
%iterative section
[rr,cc]=size(im);
[xm,ym]=TrackXY_by_COM_2Dmoment(im); 
xnw=xm;
ynw=ym;

errx=zeros(QI.iterations,1);
prequit=0;
for ii=1:QI.iterations   
     if ~prequit %Note that this may stop the loop prematurely 
        xol=xnw;
        yol=ynw;
        Xsamplinggrid=QI.X0samplinggrid+xnw;
        Ysamplinggrid=QI.Y0samplinggrid+ynw;
        allprofiles=(interp2(im,Xsamplinggrid,Ysamplinggrid,'NaN'));
        
        [aa,rara]=size(Xsamplinggrid);
        spokesnoperquad=round(aa/4);
        Qiprofs=zeros(4,rara);
        Qiprofs(1,:)=nanmean(allprofiles(1:spokesnoperquad,:));  %East
        Qiprofs(2,:)=nanmean(allprofiles(spokesnoperquad+1:2*spokesnoperquad,:));  %North
        Qiprofs(3,:)=nanmean(allprofiles(2*spokesnoperquad+1:3*spokesnoperquad,:));  %West
        Qiprofs(4,:)=nanmean(allprofiles(3*spokesnoperquad+1:4*spokesnoperquad,:));  %South

        QiHor=[fliplr(Qiprofs(3,:)) Qiprofs(1,:)];
        QiVer=[fliplr(Qiprofs(4,:)) Qiprofs(2,:)];


        %Get-image centered position, corrected for oversampling and off-center
        %sampling
        fudgefactor=(pi/2);
        xnw=-((length(QiHor)/2-SymCenter(QiHor))+0.5)/QI.radialoversampling/fudgefactor+xnw;
        ynw=-((length(QiVer)/2-SymCenter(QiVer))+0.5)/QI.radialoversampling/fudgefactor+ynw;
        errx(ii)=((xnw-xol)^2+(ynw-yol)^2);
        if (isnan(xnw)|isnan(ynw)) 
            prequit=1;
        end  
     else
         prequit=1;
         xnw=xol;  %fetch the old values
         ynw=yol;
     end  
     QI.xpos=xnw;
     QI.ypos=ynw;
end

if  sho | (nargin <3)    
     [aa,rara]=size(Xsamplinggrid);
     spokesnoperquad=round(aa/4);     
     figure(3);
     subplot(2,2,1);
         hold off
         pcolor(im); colormap bone; shading flat; hold on;
         plot(Xsamplinggrid(1:spokesnoperquad,:)'+0.5,Ysamplinggrid(1:spokesnoperquad,:)'+0.5,'r-');
         plot(Xsamplinggrid(spokesnoperquad+1:2*spokesnoperquad,:)'+0.5,Ysamplinggrid(spokesnoperquad+1:2*spokesnoperquad,:)'+0.5,'k-');
         plot(Xsamplinggrid(2*spokesnoperquad+1:3*spokesnoperquad,:)'+0.5,Ysamplinggrid(2*spokesnoperquad+1:3*spokesnoperquad,:)'+0.5,'r-');
         plot(Xsamplinggrid(3*spokesnoperquad+1:4*spokesnoperquad,:)'+0.5,Ysamplinggrid(3*spokesnoperquad+1:4*spokesnoperquad,:)'+0.5,'k-');
         plot(xm+0.5,ym+0.5,'wx', 'Markersize', 15);
         plot(xol+0.5,yol+0.5,'w+', 'Markersize', 15);
         plot(xnw+0.5,ynw+0.5,'wsq', 'Markersize', 15);     
         title('Image&sampling grid') ;
     
     subplot(2,2,2);
         hold off
         plot(Qiprofs'); hold on;
         title('Quadrant profiles');
         legend('East', 'North' , 'West', 'South');
         xlabel('position, sampling units')
         ylabel('singal, a.u')
     
     subplot(2,2,3);
         hold off
         plot(QiHor, 'r-'); hold on;
         plot(QiVer);
         title('concatenated')
          legend('Hor', 'Ver');
         xlabel('position, sampling units')
         ylabel('singal, a.u')
          
     subplot(2,2,4);
         hold off
         plot(errx,'k-o'); hold on;
         title('diff error per iteration');
         xlabel('iteration')
         ylabel('diff. error, pixels')
      
end


%Cell edge detection--------------------------------
[elliptoid,CellEdges]=Get_CellEdges(QI,allprofiles,im);

[Elliptoid_Corrected,CellEdges_Corrected, Angles_corrected,Flags]=Correct_CellEdges(QI,elliptoid,CellEdges,im);

Elliptoid_Corrected=JKD1_PRF_smooth(Elliptoid_Corrected,4);  %smooth
Elliptoid_Corrected=(Elliptoid_Corrected-nanmax(Elliptoid_Corrected))+1;  %  max is one
CellEdges_Corrected=JKD1_PRF_smooth(CellEdges_Corrected,4);  %smooth

CellEdgesVsangle=[Angles_corrected CellEdges_Corrected Elliptoid_Corrected Flags];
QI.CellEdgesVsangle=CellEdgesVsangle;

QI.PolarContourEdge_Angles=Angles_corrected;
QI.PolarContourEdge_Radii=CellEdges_Corrected;
QI.PolarContourEdge_Elliptoid=Elliptoid_Corrected;
QI.PolarContourEdge_EdgeCorrectionFlags=Flags;
QI.CartesianContourEdge_X=QI.xpos+0.5+CellEdges_Corrected.*cos(Angles_corrected);
QI.CartesianContourEdge_Y=QI.ypos+0.5+CellEdges_Corrected.*sin(Angles_corrected);



%final plot
figure(7);
dum=pcolor(im), colormap bone; hold on;    
 title('ori image plus fitted elliptoids'); xlabel('xposition, pixels'); ylabel('yposition, pixels'); axis tight; hold on;
 [rr,cc]=size(im);
 for kk=1:15
    plot(QI.xpos+0.5+kk/15*rr/2*Elliptoid_Corrected.*cos(Angles_corrected),...
         QI.ypos+0.5+kk/15*rr/2*Elliptoid_Corrected.*sin(Angles_corrected), 'w-'); hold on;
 end 
 plot(QI.xpos+0.5+CellEdges_Corrected.*cos(Angles_corrected),...
      QI.ypos+0.5+CellEdges_Corrected.*sin(Angles_corrected), 'r-', 'Linewidth', 4); hold on;




function [Elliptoid_Corrected,CellEdges_Corrected,Angles_Corrected,Flags]=Correct_CellEdges(QI,elliptoid,CellEdges,im);
%Manual correction of section
   figure(4);
    [rr,cc]=size(im);
    cor_x=500/cc;
    cor_y=500/rr;
        dum=P_Color(im,500,500, 'bone'); hold on;    
         title('ori image plus fitted elliptoids'); xlabel('xposition, pixels'); ylabel('yposition, pixels'); axis tight; hold on;
         for kk=1:15
            plot(QI.xpos*cor_x+0.5+kk/15*rr/2*elliptoid.*cos(QI.angles)*cor_x,QI.ypos*cor_y+0.5+kk/15*rr/2*elliptoid.*sin(QI.angles)*cor_y, 'w-'); hold on;
         end 
         plot(QI.xpos*cor_x+0.5+CellEdges.*cos(QI.angles)*cor_x,QI.ypos*cor_y+0.5+CellEdges.*sin(QI.angles)*cor_y, 'r-', 'Linewidth', 4); hold on;
         
          
         but=1;
         xuc=[];
         yuc=[];
         while but==1;
             [xi,yi,but]=ginput(1);
             if but==1
                xuc=[xuc; xi]; 
                yuc=[yuc; yi]; 
             end
         end
         x=xuc/cor_x; 
         y=yuc/cor_y;
     
     Flags=ones(length(CellEdges),1);
     if ~isempty(xuc)
            plot(xuc,yuc,'bo-');
            %Calculate baCK TO REGULAR COORDINATES 
            EdgContactAngles=atan2((y-QI.ypos),(x-QI.xpos));
            EdgContactRadii=((x-QI.xpos).^2+ (y-QI.ypos).^2).^0.5;

          %some cleaning up
          EdgContactAngles=mod(EdgContactAngles,2*pi);
          [EdgContactAngles,ix]=sort(EdgContactAngles);
          EdgContactRadii=EdgContactRadii(ix);

          angles=mod(QI.angles,2*pi);
          [angles,ix]=sort(angles);
          CellEdges=CellEdges(ix);

          %replacement 
          sel=find(angles<=EdgContactAngles(1) | ...
                    (angles>=EdgContactAngles(end)));
          Flags(sel)=0;
      Angles2Keep1=[angles(sel); EdgContactAngles];
      [Angles2Keep,ix]=sort(Angles2Keep1);
      Radii2Keep1= [CellEdges(sel); EdgContactRadii]; 
      Radii2Keep=Radii2Keep1(ix);
     else
          angles=mod(QI.angles,2*pi);
          [angles,ix]=sort(angles);
          Angles2Keep=angles;
          Radii2Keep=CellEdges(ix);
     end
figure(5)
plot(Angles2Keep,Radii2Keep,'ro'); hold on;
if ~isempty(xuc), plot(EdgContactAngles,EdgContactRadii,'bo-'); hold on; end
plot(Angles2Keep,Radii2Keep,'k-'); hold on;
xlabel('angle, radians');
ylabel('Cell radius');
legend('automatic', 'corrected', 'to keep');
dum=1;




CellEdges_Corrected=interp1(Angles2Keep,Radii2Keep, angles,'linear')
Elliptoid_Corrected=CellEdges_Corrected/max(CellEdges_Corrected);  %normalized
Angles_Corrected=angles;  




function [elliptoid,CellEdges]=Get_CellEdges(QI,allprofiles,im);
[rr,cc]=size(allprofiles);
CellEdges=zeros(rr,1);
[~,ix]=max(std(allprofiles'));  %find most prominent edge
refprf=allprofiles(ix,:);
%refprf=nanmean(allprofiles);
%[AvCellEdge,~, ~, ~,~,~, ~,~]=Splitfast(refprf); 
[~,AvCellEdge]=min(refprf);

for ii=1:rr
    prf=allprofiles(ii,:);
    mp=nanmean(prf);
    sel=find(isnan(prf)); prf(sel)=mp;  %padding nans
    fw=prf/nanstd(prf);   fw=fw-nanmean(fw);            
    rv=refprf/nanstd(refprf);  rv=rv-nanmean(rv);      
    d=real(ifft(fft(fw).*conj(fft(rv))));
    ld=ceil(length(d)/2);
    d=[d(ld+1:length(d)) d(1:ld)]';   %swap first and second half 
    %plot(d); pause(0.5);
    [~,x]=max(d);
    shft=(x-ld);  %this is now the  shift of the profile relative to the reference profile
    %x=(subpix_step(d)); 
    CellEdges(ii)=(shft+AvCellEdge);    
end

elliptoid=CellEdges/max(CellEdges);  %normalized
CellEdges=CellEdges/QI.radialoversampling;



data=[QI.angles CellEdges]; 
CellEdges_Expanded=Expand_elliptoidrange(data);




function x=subpix_step(d);
    %this function performs a subpixel step by parabolic fitting
    hf=3; ld=length(d);  xs=[1:1:ld]';   [val,x]=min(d);  %uneven hf
    lo=max([x-hf 1]); hi=min([x+hf ld]);  %cropping
    ys=d(lo:hi); xs=xs(lo:hi);
    prms=polyfit(xs,ys,2); x=-prms(2)/(2*prms(1));



function outdata=JKD1_PRF_smooth(data,span)
%running window 'span' smoothing average; 0=no smoothing.
%perioc boundaries
%JacobKers 2013

if span>0
    hs=ceil(span/2);
    le=length(data);
    outdata=zeros(le,1);
    midpt=ceil(le/2);
    data2=[data(midpt:end) ; data ; data(1:midpt)];
    lex=length(data(midpt:end));
    for i=1:le
    idx=i+lex;
    outdata(i)=mean(data2(idx-hs:idx+hs));
    end
else
    outdata=data; %no smoothing
end

function data_expanded=Expand_elliptoidrange(data);
    %Expand for later interpolation
    [rr,~]=size(data);
    halfidx=ceil(rr/2);
    half1=data(1:halfidx-1,:);
    half2=data(halfidx:end,:);
    half1(:,1)=half1(:,1)+2*pi;
    half2(:,1)=half2(:,1)-2*pi;
    data_expanded=[half2 ; data ; half1];



           
 
       



