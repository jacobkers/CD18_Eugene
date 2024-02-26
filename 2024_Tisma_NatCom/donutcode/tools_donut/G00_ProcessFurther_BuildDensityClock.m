function G00_ProcessFurther_BuildDensityClock(DistAxisPerc,AvContentDist,initval);;
%This function builds a visual representation of (a number of) chromosome
%maps

    close all;
    
    action.BuildSecondaryRings=1;
    
    HalfSize=50;
    InnerCircle=HalfSize/2;
    %CircleRange=InnerCircle*0.38/1.600;
    CircleRange=HalfSize/8;
    PlotTheWidth=0;
    RotAngle=0/180*pi; %this rotates the distance&colormap relative to the genome map
  
    Canvas=zeros(2*HalfSize+1)+max(AvContentDist);
    [XX,YY]=meshgrid(-HalfSize:HalfSize,-HalfSize:HalfSize);
    RR=((XX.^2)+(YY.^2)).^0.5;
    Arg=atan2(YY,XX);
    
    Content_Norm=CircleRange*AvContentDist/range(AvContentDist);
    
    %build a mapped axis
    GG=cumsum(AvContentDist); 
    GenAxisPerc=100*GG/max(GG);
     
    %Set inner and outer Radii
    RGen=0*Content_Norm+0.8*InnerCircle;
    RIn=0*Content_Norm+InnerCircle;
    
    if PlotTheWidth
        ROut=InnerCircle+Content_Norm;
    else
        ROut=InnerCircle+0*Content_Norm+CircleRange;
    end
    Angle0=(100-initval.Rfppos)/100*2*pi+0.5*pi;
    
   [OriX, OriY, DifX, DifY, RfpX,RfpY,CfpX,CfpY,TicksX,TicksY]=Build_Markers(RGen(1),initval,HalfSize);
    
        
    
    [GenSpokesX, GenSpokesY,GenX,GenY]=MapDistanceToGenome(Angle0,RotAngle,DistAxisPerc,GenAxisPerc,...
                                                            HalfSize,RIn,RGen);        
    %build ring 1
    [Canvas,InX,InY,OutX,OutY]=ColorHeatMapRing(Canvas, Angle0+RotAngle,DistAxisPerc,...   %inner ring
      AvContentDist,HalfSize,RIn,ROut);  
    
  
  if action.BuildSecondaryRings
    inpth='D:\jkerssemakers\My Documents\BN CD Recent\BN_CD15 FabaiXuanCells\AnalysisResults_WideField\RingData\';
    HU=load(strcat(inpth,'AllCells_HU.mat'),'AvContentDist');
    AvContentDist_HU=HU.AvContentDist;
    Fis=load(strcat(inpth,'AllCells_Fis.mat'),'AvContentDist');
    AvContentDist_Fis=Fis.AvContentDist;
    HNS=load(strcat(inpth,'AllCells_HNS.mat'),'AvContentDist');
    AvContentDist_HNS=HNS.AvContentDist;

  %Build Ring 2
     RIn2=ROut+CircleRange/5;     
     ROut2=RIn2+CircleRange/2;
    [Canvas,InX2,InY2,OutX2,OutY2]=ColorHeatMapRing(Canvas, Angle0+RotAngle,DistAxisPerc,...
                                                AvContentDist_HU,HalfSize,RIn2,ROut2);
                                            
     %Build Ring 3
     RIn3=ROut2+CircleRange/5;     
     ROut3=RIn3+CircleRange/2;
    [Canvas,InX2,InY2,OutX2,OutY2]=ColorHeatMapRing(Canvas, Angle0+RotAngle,DistAxisPerc,...
                                                AvContentDist_Fis,HalfSize,RIn3,ROut3);
                                            
      %Build Ring 4
     RIn4=ROut3+CircleRange/5;     
     ROut4=RIn4+CircleRange/2;
    [Canvas,InX2,InY2,OutX2,OutY2]=ColorHeatMapRing(Canvas, Angle0+RotAngle,DistAxisPerc,...
                                                AvContentDist_HNS,HalfSize,RIn4,ROut4);
 
  
  
  end
  
  %plot the result
    [rr,cc]=size(Canvas);
    P_Color(Canvas,cc,rr,'hot'); colormap hot; shading flat; axis equal; hold on; 
 
    plot(InX,InY,'k', 'LineWidth',1);
    plot(OutX,OutY,'k', 'LineWidth',2);
    plot(GenSpokesX,GenSpokesY,'k-', 'LineWidth',1);
    
    plot(GenX,GenY,'k', 'LineWidth',2); hold on;
    plot(RfpX,RfpY,'ro','MarkerFacecolor','r', 'MarkerSize',5); %,CfpX,CfpY,TicksX,TicksY
    plot(CfpX,CfpY,'co','MarkerFacecolor','c', 'MarkerSize',5);
    
    plot(OriX,OriY,'ko', 'MarkerFacecolor','k', 'MarkerSize',5); hold on;
    plot(DifX,DifY,'ko', 'MarkerFacecolor','k', 'MarkerSize',5); hold on;
    
    
    
    
    plot(TicksX,TicksY,'ko','MarkerFacecolor','k', 'MarkerSize',5);
    
    saveas(gcf,strcat(inpth,'Rings.tif'),'tif');  
    
    function nearpoints=get_near_points(x0,y0, r0, points)
    %this function returns  coordinates within a distance r0 from a point x0,y0 in a 'peaks'(x,y,I,) database
    %JacobKers 2012
    rall=((points(:,1)-x0).^2+(points(:,2)-y0).^2).^0.5;
    %sel=find(rall<r0);                
    nearpoints=points(rall<r0,:);
    
 function [GenSpokesX, GenSpokesY,GenX,GenY]=MapDistanceToGenome(Angle0,RotAngle,DistAxisPerc,GenAxisPerc,...
                                                           HalfSize,RIn,RGen);
    %Build spokes linking spatial distance and genomic location
    skips=5;
    AngleRun=Angle0+RotAngle+2*pi*(-DistAxisPerc/100);
    InX=HalfSize+1+RIn.*cos(AngleRun);       %inner border ring #1
    InY=HalfSize+1+RIn.*sin(AngleRun);
    AngleRunGen=Angle0+2*pi*(-GenAxisPerc/100);
    GenX=HalfSize+1+RGen.*cos(AngleRunGen);  %genome map #1
    GenY=HalfSize+1+RGen.*sin(AngleRunGen);
    GenSpokesX=[GenX(1:skips:end); InX(1:skips:end)]; 
    GenSpokesY=[GenY(1:skips:end); InY(1:skips:end)];

    
 function [Canvas,InX,InY,OutX,OutY]=ColorHeatMapRing(Canvas, Angle0,DistAxisPerc,...
                                                AvContentDist,HalfSize,RIn,ROut);
            %Build a heat map ring 
    LL=length(DistAxisPerc);        
    %Determine Cartesian positions of the heat map circle
    AngleRun=Angle0+2*pi*(-DistAxisPerc/100);
    InX=HalfSize+1+RIn.*cos(AngleRun);       %inner border ring #1
    InY=HalfSize+1+RIn.*sin(AngleRun);
    
    OutX=HalfSize+1+ROut.*cos(AngleRun);    %outer border ring #1
    OutY=HalfSize+1+ROut.*sin(AngleRun);
       
    %Build a radial 'spoke' grid in the circle area
    ClockGridX=zeros(LL,10); ClockGridY=zeros(LL,10);    
    for jj=1:LL
        ClockGridX(jj,:)=linspace(InX(jj),OutX(jj),10);
        ClockGridY(jj,:)=linspace(InY(jj),OutY(jj),10);
    end
       
    %Color the canvas
    [XX,YY]=meshgrid(-HalfSize:HalfSize,-HalfSize:HalfSize);
    points=[XX(:) YY(:)]+HalfSize;
    for jj=1:LL
        spokeX=ClockGridX(jj,:); spokeY=ClockGridY(jj,:);
        localRIn=RIn(jj); localROut=ROut(jj);
        for ii=1:10
            spokeXi=spokeX(ii); spokeYi=spokeY(ii);
            nearpoints=get_near_points(spokeXi,spokeYi, 3, points); %nearby points
            LocalR=((nearpoints(:,1)-HalfSize).^2+(nearpoints(:,2)-HalfSize).^2).^0.5;
            sel=find((LocalR>localRIn)&(LocalR<localROut)); %crop
            nearpoints=nearpoints(sel,:);            
            [LN,~]=size(nearpoints);
            for kk=1:LN
            Canvas(nearpoints(kk,2),nearpoints(kk,1))=AvContentDist(jj);
            end
        end
    end
    
    function [OriX, OriY, DifX, DifY, RfpX,RfpY,CfpX,CfpY,TicksX,TicksY]=...
            Build_Markers(Radius,initval,HalfSize);
        %Build some genomic markers
        AngleRfp=(100-initval.Rfppos)/100*2*pi+0.5*pi;
        AngleCfp=(100-initval.Cfppos)/100*2*pi+0.5*pi;
        AngleOri=(100-initval.Oripos)/100*2*pi+0.5*pi;
        AngleDif=(100-initval.Difpos)/100*2*pi+0.5*pi;
        
        
        
        AngleTicks=linspace(0,2*pi,5);
        
        RfpX=HalfSize+1+Radius*cos(AngleRfp);       %inner border ring #1
        RfpY=HalfSize+1+Radius*sin(AngleRfp);
        
        CfpX=HalfSize+1+Radius*cos(AngleCfp);       %inner border ring #1
        CfpY=HalfSize+1+Radius*sin(AngleCfp);
        
        OriX=HalfSize+1+Radius*cos(AngleOri);       %inner border ring #1
        OriY=HalfSize+1+Radius*sin(AngleOri);
        
        DifX=HalfSize+1+Radius*cos(AngleDif);       %inner border ring #1
        DifY=HalfSize+1+Radius*sin(AngleDif);
        
        TicksX=HalfSize+1+Radius*cos(AngleTicks);
        TicksY=HalfSize+1+Radius*sin(AngleTicks);
        