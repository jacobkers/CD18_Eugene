    
function [xnw,ynw,allprofiles]=TrackXY_by_QI(im,QI,sho);
%This function prepares a sub-pixel XY fit by making 4 profiles in QI
%style
  if nargin <3
      sho=1;
    close all
    test.PicSize=50; 
    x0=test.PicSize/2+17; 
    y0=test.PicSize/2; 
    test.PatternRingRadius=test.PicSize/8;
    test.PatternRingWidth=test.PicSize/30;
    im=MakeHighResRing(x0,y0,test); 
    
    QI.radialoversampling=2;
    QI.angularoversampling=0.7;
    QI.minradius=0;
    QI.maxradius=test.PicSize/4;
    QI.iterations=10;
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
        
end

if sho %nargin <3  %Plotting Menu
    
     [aa,rara]=size(Xsamplinggrid);
     spokesnoperquad=round(aa/4);
          
     subplot(2,2,1);
     hold off
     pcolor(im); colormap bone; shading flat; hold on;
     skipspokes=10;
     plot(Xsamplinggrid(1:skipspokes:spokesnoperquad,:)'+0.5,Ysamplinggrid(1:skipspokes:spokesnoperquad,:)'+0.5,'bo-', 'MarkerSize', 1);
     plot(Xsamplinggrid(spokesnoperquad+1:skipspokes:2*spokesnoperquad,:)'+0.5,Ysamplinggrid(spokesnoperquad+1:skipspokes:2*spokesnoperquad,:)'+0.5,'b-o', 'MarkerSize', 2);
     plot(Xsamplinggrid(2*spokesnoperquad+1:skipspokes:3*spokesnoperquad,:)'+0.5,Ysamplinggrid(2*spokesnoperquad+1:skipspokes:3*spokesnoperquad,:)'+0.5,'r-o', 'MarkerSize', 2);
     plot(Xsamplinggrid(3*spokesnoperquad+1:skipspokes:4*spokesnoperquad,:)'+0.5,Ysamplinggrid(3*spokesnoperquad+1:skipspokes:4*spokesnoperquad,:)'+0.5,'b-o', 'MarkerSize', 2);
     plot(xm+0.5,ym+0.5,'wx', 'Markersize', 15);
     plot(xol+0.5,yol+0.5,'w+', 'Markersize', 15);
     plot(xnw+0.5,ynw+0.5,'wsq', 'Markersize', 15);
     
     title('Image&sampling grid') ;
          
     subplot(2,2,2);
     hold off
     plot(Qiprofs'); hold on;
     title('profiles');
     legend('East', 'North' , 'West', 'South');
     
     subplot(2,2,3);
     hold off
     plot(QiHor, 'r-'); hold on;
     plot(QiVer);
     title('concatenated')
      legend('Hor', 'Ver');
          
     subplot(2,2,4);
     hold off
     plot(errx,'k-o'); hold on;
     [~]=ginput(1);      close(gcf)
     pause(0.05);

end

    
 
       



