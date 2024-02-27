function A005_Build_MovieList(batchrunindex)
%JWJK_A:-------------------------------------------------------------------
%Build movielist

%Description: identifies sets of frames per movie based on the filenames

%Reference: CD lab, project Sandro, written by Jacob Kers 2018-20
%:JWJK_A-------------------------------------------------------------------

close all;
if nargin<1,batchrunindex=5;end
    
initval=A000_Repli_Init(batchrunindex);

%First, make a list of cells via the first frame
[CellNames,timedigits]=Build_Cell_list(initval);

NCells=length(CellNames);
disp(strcat(initval.expi,'  Counted',num2str(NCells),'movies'));

for cc1=1:NCells  %for all cells:            
    CellName=char(CellNames(cc1));
    CellName=CellName(1:8);  %cell name
    CellFrames=cell(0);    
    sel=strfind(initval.Cell_Labels,CellName);
    for ii=1:length(sel);
        if ~isempty(sel{ii})
            ThisCellFrameName=initval.Cell_Labels{ii};
            CellFrames=[CellFrames ThisCellFrameName];            
            FrameTimes(ii)=str2num(char(ThisCellFrameName(timedigits)));
            %CellName1=strcat('ResultsOfCell',ThisCellFrameName,'.mat'); 
%            MatFilePath=strcat(initval.resultpath,'ResultsPerCellMatlab',initval.DirSep);
%            load(strcat(MatFilePath,CellName1));        
        end
    end
    MovieList(cc1).CellNames=CellName;  %name of cell
    MovieList(cc1).CellFrames=CellFrames; %names of cell frames
    
    if initval.onlyfirstmovieframe,
        MovieList(cc1).FrameTimes=FrameTimes(1);
    else
        MovieList(cc1).FrameTimes=FrameTimes;
    end
    
end


MatFilePath=strcat(initval.resultpath,'ResultsPerCellMatlab',initval.DirSep);
ResultName=strcat(MatFilePath,strcat('__MovieList.mat'));
save(ResultName, 'MovieList');

function [CellNames,timedigits]=Build_Cell_list(initval);
CellNames=cell(0);
LC=length(initval.Cell_Labels);
timedigits=[];
for ii=1:LC
    ThisCellName=initval.Cell_Labels{ii};
    istwo=strfind(ThisCellName,'t01');
    isone=strfind(ThisCellName,'t1');
    if ~isempty(istwo),
        CellNames=[CellNames initval.Cell_Labels{ii}];
        timedigits=10:11;
    end;
    if ~isempty(isone),
        CellNames=[CellNames initval.Cell_Labels{ii}];
        timedigits=10;
    end;  
end
dum=1;

    
    
    
    
