function CellEdgedistFromChomosomeCenter=F005_Relate_CromosomeToCellEdge(initval, savit)
% 'Use this section for a Quicksheet'
%------------------------------------------------------------------------
%This function : %links chromosome position to cell wall position; 
%we use two definitions: 
%1) the radial distancefrom the chromosome center (COM-based on contour) 
    %to the first cell wall crossing
%2) %the nearest xy cell wall position for a chromosome xy position

%The results are added as fileds to the 'Chromosome' data structure as:
        % RadialCellEdge: [252x1 double]
        % CellEdgeNearestDist: [1x252 double]
        % CellEdgeNearestIndex: [1x252 double]
%------------------------------------------------% Jacob Kerssemakers 2016
% 'End of Quicksheet section'
load(initval.PathResultName);    
[rr,cc,dd]=size(yfp_Stack);
 h=figure;
 
%% Analysis section; data is stored
    close all;        
 for jj=1:dd  
       Chromosome=AllFrameData(jj).Chromosome; 
       BrightField=AllFrameData(jj).BrightField;
       
       NoGaps=find((1.0-Chromosome.Gaps)==1);
       CHangle=Chromosome.AnnularAxis;
       CHx0=Chromosome.xCOM;
       CHy0=Chromosome.yCOM;
       BX=2*BrightField.CartesianContourEdge_X;
       BY=2*BrightField.CartesianContourEdge_Y;
       CX=Chromosome.CartesianContourEdge_X;
       CY=Chromosome.CartesianContourEdge_Y;
       
       
       %% 1 get radial cell edge distance to chromosome center, sorted 0-2*pi
       BangleC1=atan2((BY-CHy0),(BX-CHx0));  %angles
       BangleC2=mod(BangleC1,2*pi);
       [BangleC,ix]=sort(BangleC2);            %sorted 0-2pi   
       BradC=((BX-CHx0).^2+(BY-CHy0).^2).^0.5; %radial distances
       BradC=BradC(ix);                        %sorted 0-2pi 
       CellEdgedistFromChomosomeCenter=interp1(BangleC, BradC,CHangle);  %radii mapped on chromosome angles
       AllFrameData(jj).Chromosome.RadialCellEdge=CellEdgedistFromChomosomeCenter;
       %% 2 get the nearest cell wall point per chromosome position
       LC=length(CX);
       CellEdgedistNearest=zeros(LC,1);
       CellEdgedistNearestIndex=zeros(LC,1);
       for cc=1:LC  %for  all cell positions
           cx=CX(cc);
           cy=CY(cc);
           [mindist,ix]=min(((BX-cx).^2+(BY-cy).^2).^0.5); %minimum radial distance
           CellEdgeNearestDist(cc)=mindist;   %minimum distance
           CellEdgeNearestIndex(cc)=ix;       %minimum point (index)
            %minimum edge
       end
       AllFrameData(jj).Chromosome.CellEdgeNearestDist=CellEdgeNearestDist;
       AllFrameData(jj).Chromosome.CellEdgeNearestIndex=CellEdgeNearestIndex;      
       %% 
       
       if 0
           close all
           plot(CX(NoGaps),CY(NoGaps),'ro-'); hold on;
           plot(BX,BY,'ko-'); hold on;
           NBX=BX(CellEdgeNearestIndex);
           NBY=BY(CellEdgeNearestIndex);
           pairs1=[CX(NoGaps)';NBX(NoGaps)']; 
           pairs2=[CY(NoGaps)';NBY(NoGaps)'];
           plot(pairs1,pairs2,'-'); hold on;           
           [~]=ginput(1);
       end
 end
    
%----------------------------------------------
if savit,save(initval.PathResultName,'AllFrameData','/append'); end


