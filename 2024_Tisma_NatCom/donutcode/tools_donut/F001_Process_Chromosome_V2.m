function F001_Process_Chromosome_V2(initval,actions)
% 'Use this section for a Quicksheet'
%-------------------------------------------------------------------------
%Sequence of actions to work through chromosome
%This function adds the following fields to the data structure 'Chromosome'
    % Chromosome = 
        % RadialMap: [252x80 double]
        % AnnularAxis: [252x1 double]
        % RadialAxis: [1x80 double]
        % MappingIndex: [1x252 double]
        % PolarContourMax: [252x1 double]
        % PolarContourEdge: [252x1 double]
        % PolarContourRim: [252x1 double]
        % PolarContourContent: [1x252 double]
        % PolarContourFWHM: [1x252 double]
        % Gaps: [1x252 logical]
        % xCOM: 102.6826
        % yCOM: 100.0530
        % CartesianContourMax_X: [252x1 double]
        % CartesianContourMax_Y: [252x1 double]
        % CartesianContourEdge_X: [252x1 double]
        % CartesianContourEdge_Y: [252x1 double]
        % picture: [201x221 double]
 %Later more fields are added:(Function F005_Relate_CromosomeToCellEdge) 
        % RadialCellEdge: [252x1 double]       
        % CellEdgeNearestDist: [1x252 double]
        % CellEdgeNearestIndex: [1x252 double]
%------------------------------------------------% Jacob Kerssemakers 2016 
% 'End of Quicksheet section'


load(initval.PathResultName);    
[rr,cc,dd]=size(yfp_Stack);
 h=figure;
 
%% Analysis section; data is stored
if actions.yfp_analyze
    close all;    
    
 for jj=1:dd
    dd-jj+1
    pic=yfp_Stack(:,:,jj);      
    dum=1;
        CellSoftMask=AllFrameData(jj).BrightField.CellMask2X;
        maskpic=CellSoftMask.*(pic- median(pic(:)));
        QI_chro=TrackXY_by_QI_Init_Chro(maskpic);  
        
        Chromosome=ChromosomeRadialSampling(maskpic,QI_chro, 0,initval);       
        Chromosome.picture=pic; 
        
        [aa,rr]=size(Chromosome.RadialMap);
        
        if 0
            figure(4);
            MappingIndex=(1:aa);
            VI=Chromosome.ValidIdx;
            cr=QI_chro.radialoversampling
            pcolor(cr*Chromosome.RadialMap); colormap bone; shading flat; hold on;
            plot(cr*Chromosome.PolarContourMax(VI), MappingIndex(VI), 'ro', 'Linewidth', 2);
            plot(cr*Chromosome.PolarContourEdge(VI), MappingIndex(VI), 'bo', 'Linewidth', 1);
            plot(cr*Chromosome.PolarContourRim(VI), MappingIndex(VI), 'bo', 'Linewidth', 1);
            title('radial map of brightfield with bf-edge overlay');
            xlabel('radial sampling steps');
            ylabel('angular sampling steps');
            [~]=ginput(1);
        end
        AllFrameData(jj).Chromosome=Chromosome; 
    end
    
%----------------------------------------------
end
save(initval.PathResultName,'AllFrameData','/append');


