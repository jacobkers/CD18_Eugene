function A010_WF_PerCell_AnalyzeCellShapeStandAlone(initval)
%JWJK_A:-------------------------------------------------------------------
%Description: obtains general cell shape data

%input: .mat database from crop code analysis.

%output:.mat database is generated and expanded at later steps.

%Reference: CD lab, project Sandro, written by Jacob Kers 2018-20
%:JWJK_A-------------------------------------------------------------------


close all;
%runtime options
actions.workfullstackdiagnosis=0;
actions.showandsavepics=1;
if nargin<1
    usr='Jacob', batchrunindex=0;
    initval=A000_Repli_Init(batchrunindex,usr);
end

disp(initval.expi);
CellImagePath=strcat(initval.pth_repli,'CellImages_All',initval.DirSep);
if isdir(CellImagePath), 
    rmdir(CellImagePath,'s');  
end
mkdir(CellImagePath);
MatFilePath=strcat(initval.pth_repli,'ResultsPerCellMatlab',initval.DirSep);
LC=length(initval.Cell_Labels);

for ii=1:LC
    cellno=char(initval.Cell_Labels{ii});    
    CellName=strcat('ResultsOfCell',cellno); 
    disp(strcat('Program:A10_experiment:',initval.expi,':',CellName,'CellshapeAnalysis..', num2str(LC-ii+1), 'cells to go'));   
    NumCellLabel=BuildNumericCellLabel(cellno);   
    
    %load and identify color channels:
    channel_stack=get_channel_stack_from_cropcode(cellno,initval,actions); 
    [cellmask, ~, ~, ~, ~]=get_channel_ID(channel_stack, initval);
     
    edge_pic = bwmorph(cellmask,'remove');      
    Cell = regionprops('struct',cellmask,...
        'Centroid','MajorAxisLength','MinorAxisLength',...
        'Area','Orientation', 'EquivDiameter');
    Cell.BW=cellmask;
    Cell.Edge=edge_pic;  
    Cell.EdgeLength=sum(1*edge_pic(:));
    Cell.PerimeterFactor=Cell.EdgeLength/(Cell.EquivDiameter*pi);
    %Make course estimate of number of cells    
    N_est=Estimate_Cellnumber(Cell,initval);    
    Cell.N_est=N_est;    
    initval.ResultName=strcat(MatFilePath,strcat('ResultsOfCell',cellno,'_Cellshape.mat'));
    save(initval.ResultName, 'Cell');
end

function pic=GetWorkpicFromStack(stack,WorkPicOption);
    stack=double(stack);
    switch WorkPicOption
        case 'MeanProject',pic=squeeze(mean(stack,3));
        case 'FocalPlane', 
        %get main plane  (assuming largely planar features)
        stcurve=squeeze(std(std(stack)));
        [~,MainPlane]=max(stcurve);
        pic=double(stack(:,:,MainPlane));  
        case 'First', pic=double(stack(:,:,1));
    end
       
function N_est=Estimate_Cellnumber(Cell,initval);
    ballparksinglecellarea=2000*2000/((initval.nmperpixel)^2); %pixels squared
    ballparksinglecellperimeter=2*pi*((ballparksinglecellarea/pi)^0.5);
    N_est_early=(Cell.Area/ballparksinglecellarea);   
    N_est=round(N_est_early)*1.3;
    dum=1;
 
       
