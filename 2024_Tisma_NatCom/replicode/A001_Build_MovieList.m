function A001_Build_MovieList
%JWJK_A:-------------------------------------------------------------------
%Build movielist

%Description: identifies sets of frames per movie based on the filenames

%Reference: CD lab, project Sandro, written by Jacob Kers 2018-20
%:JWJK_A-------------------------------------------------------------------
%find movieframes
close all;

batchrunindex=102;
initval=A000__WF_Get_JacobPathsandExperimentsReplicode(batchrunindex);

CellNames=cell(0);
%First, make a list of cells via the first frame
sel=strfind(initval.Cell_Labels,'t01');
for ii=1:length(sel);
    if ~isempty(sel{ii})
        CellNames=[CellNames initval.Cell_Labels{ii}];
    end
end
NCells=length(CellNames);


for cc1=1:NCells  %for all cells:        
    disp(strcat('Correlation Analysis..',num2str(NCells-cc1+1), 'cells to go'));
    CellName=char(CellNames(cc1));
    CellName=CellName(1:8);  %cell name
    CellFrames=cell(0);    
    sel=strfind(initval.Cell_Labels,CellName);
    for ii=1:length(sel);
        if ~isempty(sel{ii})
            ThisCellFrameName=initval.Cell_Labels{ii};
            CellFrames=[CellFrames ThisCellFrameName];
            FrameTimes(ii)=str2num(char(ThisCellFrameName(10:11)));
            %CellName1=strcat('ResultsOfCell',ThisCellFrameName,'.mat'); 
%            MatFilePath=strcat(initval.resultpath,'ResultsPerCellMatlab',initval.DirSep);
%            load(strcat(MatFilePath,CellName1));        
        end
    end
    MovieList(cc1).CellNames=CellName;  %name of cell
    MovieList(cc1).CellFrames=CellName; %names of cell frames
    MovieList(cc1).FrameTimes=FrameTimes;
end
MatFilePath=strcat(initval.resultpath,'ResultsPerCellMatlab',initval.DirSep);
ResultName=strcat(MatFilePath,strcat('__MovieList.mat'));
save(ResultName, 'MovieList');



    
    
    
    
