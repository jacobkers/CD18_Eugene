function F002_Process_BrightField(initval,actions);
% 'Use this section for a Quicksheet'
%-------------------------------------------------------------------------
%This function tracks the brightfield image edges and stores for later use.
%Option for manual contour corrections included
% The field of 'BrightField' = 
%     CartesianContourEdge_X: [252x1 double]
%     CartesianContourEdge_Y: [252x1 double]
%                  QICenterX: 51.9848
%                  QICenterY: 53.5217
%                   CellMask: [101x111 double]
%                 CellMask2X: [201x221 double]
%------------------------------------------------% Jacob Kerssemakers 2016 
% 'End of Quicksheet section'

if actions.bf_analyze    
    load(initval.PathResultName);    
    [rr,cc,dd]=size(bf_Stack);
    driftx0=zeros(dd,1);
    drifty0=zeros(dd,1);
    QI=TrackXY_by_QI_Init(bf_pic);
    if 1
        [Pre_Mask,OutsideArea]=PreMaskit(bf_pic,initval);
        [x0,y0,~]=TrackXY_by_QI(Pre_Mask.*bf_pic+OutsideArea,QI,0);
        BrightFieldData=CheckEllipticity_by_QI(Pre_Mask.*bf_pic+OutsideArea,QI,1);
    else
        [x0,y0,~]=TrackXY_by_QI(bf_pic,QI,0);
        BrightFieldData=CheckEllipticity_by_QI(bf_pic,QI,1);
    end
    
    
    for jj=1:dd
        dd-jj+1
        pic=bf_Stack(:,:,jj).*Pre_Mask+OutsideArea;
        %find the center of the pattern
        [xnw,ynw,~]=TrackXY_by_QI(pic,BrightFieldData,0); 
        drift.x0(jj)=xnw;
        drift.y0(jj)=ynw;        
    end 
    
     if initval.skipdriftcorrection %kill drift analysis   
        drift.x0=0*drift.x0+nanmedian(drift.x0);
        drift.y0=0*drift.y0+nanmedian(drift.y0);
     end
     
      for jj=1:dd   
         xnw=drift.x0(jj);
         ynw=drift.y0(jj);
         pic=bf_Stack(:,:,jj).*Pre_Mask+OutsideArea;
        %add averaged contour to new position
         CX=BrightFieldData.CartesianContourEdge_X-x0+xnw;
         CY=BrightFieldData.CartesianContourEdge_Y-y0+ynw;
         AllFrameData(jj).BrightField.CartesianContourEdge_X=CX;
         AllFrameData(jj).BrightField.CartesianContourEdge_Y=CY;
         AllFrameData(jj).BrightField.QICenterX=xnw;
         AllFrameData(jj).BrightField.QICenterY=ynw; 
         [CellMask,CellMask2X]=BuildCellMasks(pic,BrightFieldData,AllFrameData(jj));
         AllFrameData(jj).BrightField.CellMask=CellMask;
         AllFrameData(jj).BrightField.CellMask2X=CellMask2X;
    if 0
        close all;
        pcolor(bf_pic); colormap bone; shading flat; hold on;
        plot(xnw,ynw,'rx', 'Markersize',16);
        [~]=ginput(1);
    end
        pause(0.5);
    end
    %position correction relative to averaged image
    BrightFieldData.driftX=drift.x0-x0;
    BrightFieldData.driftY=drift.y0-y0;   
    save(initval.PathResultName,'BrightFieldData','AllFrameData','/append');
    
end

function [CellMask,CellMask2x]=BuildCellMasks(pic,BrightFieldData,ThisFrameData);
%build a mask based on the cell rim; make double-pixel version for use in
%SIM images - JK16

decay=2;
[rr,cc]=size(pic);
[XX,YY]=meshgrid(1:cc,1:rr);
xc=ThisFrameData.BrightField.QICenterX;
yc=ThisFrameData.BrightField.QICenterY;

localradii=((XX-xc).^2+(YY-yc).^2).^0.5;
localanni=atan2(YY-yc,XX-xc);
CellMask=ones(rr,cc);
CellMask2x=ones(2*rr,2*cc);
edgeangles=BrightFieldData.PolarContourEdge_Angles;
edgeradii=BrightFieldData.PolarContourEdge_Radii;
for ii=1:rr
    for jj=1:cc
        localangle=localanni(ii,jj);
        localradius=localradii(ii,jj);
        [~,ix]=min(abs(edgeangles-localangle));
        localedgeradius=edgeradii(ix);
        
        if localradius>localedgeradius  %outside the rim!
            rimdist=localradius-localedgeradius;
            maskval=exp(-((rimdist/decay)^2));
            CellMask(ii,jj)=maskval;
            CellMask2x(2*(ii-1)+1:2*(ii-1)+2,2*(jj-1)+1:2*(jj-1)+2)=maskval;
        end
    end
end
CellMask2x=CellMask2x(1:end-1,1:end-1);

function [PreMask,OutsideArea]=PreMaskit(prepic,initval);
%This function pre-masks the image
close all;
[rr,cc]=size(prepic);
OutsideArea=0*prepic
PreMask=ones(rr,cc);


cor_x=500/cc;
cor_y=500/rr;


dum=P_Color(prepic,500,500, 'bone'); hold on;    
title('ori image plus fitted elliptoids'); xlabel('xposition, pixels'); ylabel('yposition, pixels'); axis tight; hold on;
but=1;
xxi=[];
yyi=[];
while but==1;
     [xi,yi,but]=ginput(1);
     if but==1
        xxi=[xxi; xi]; 
        yyi=[yyi; yi]; 
     end
end
xxi=xxi/cor_x; 
yyi=yyi/cor_y;

lox=min(xxi); hix=max(xxi);
loy=min(yyi); hiy=max(yyi);

pic=prepic;

[XX,YY]=meshgrid(1:cc,1:rr);

BFmedlevel=median(prepic(:));
sel=find((XX<lox)|(XX>hix)|(YY<loy)|(YY>hiy));
OutsideArea(sel)=BFmedlevel;
PreMask(sel)=0;


if 0
    close all;
    pcolor(OutsideArea+PreMask.*pic); colormap bone; shading flat
    [~]=ginput(1);
end
dum=1;





