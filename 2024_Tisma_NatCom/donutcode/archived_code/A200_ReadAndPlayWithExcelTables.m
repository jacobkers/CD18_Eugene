function A200_ReadAndPlayWithExcelTables
close all;
batchrunindex=6;
initval=A000__WF_Get_JacobPathsandExperiments(batchrunindex);

[numdat,txtdat]=xlsread(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A050_ClusterColumnReport.xlsx'));
[numdat2,txtdat2]=xlsread(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A050_Cell_ClusterReport.xlsx'));

%the 'numdat' rows are shifted one, column indexes are identical
%for re-writing, use first numdat row to detrmine text or numeric input
%(for text, numeric is empty)
HeadersIN=txtdat(1,:); 
ValuesIN_num=numdat;

HeadersINPerCell_1=txtdat2(1,:); 
HeadersINPerCell_2=txtdat2(2,:); 
ValuesINPerCell=numdat2(3:end,:);

 
%% somespecific column actions
CellIndex=ValuesIN_num(1:end,find(strcmp(HeadersIN,'index')));
Contents= initval.genomelength*ValuesIN_num(1:end,find(strcmp(HeadersIN,'content')));
Areas=(initval.nmperpixel)^2*ValuesIN_num(1:end,find(strcmp(HeadersIN,'area')));
% RadGyr1=initval.nmperpixel*ValuesIN_num(1:end,find(strcmp(HeadersIN,'diameter_2Rg')));
% RadGyr2=initval.nmperpixel*ValuesIN_num(1:end,find(strcmp(HeadersIN,'diameter_ContourBW_2Rg')));
RadGyr3=initval.nmperpixel*ValuesIN_num(1:end,find(strcmp(HeadersIN,'diameter_ContourBW_Equiv')));
ClusterNoPerCell=ValuesINPerCell(1:end,find(strcmp(HeadersINPerCell_2,'cluster number')));


MxA=nanmax(Areas);
Axz=linspace(0,MxA,10);
HistArea=hist(Areas,Axz);

MxC=nanmax(Contents);
AxzC=linspace(0,MxC,10);
HistContent=hist(Contents,AxzC);

MxN=nanmax(ClusterNoPerCell);
AxzN=0:1:MxN;
HistClusterNo=hist(ClusterNoPerCell,AxzN);

subplot(1,3,1);
bar(AxzN,HistClusterNo,'r');
xlabel('nluster number per cell)')
ylabel('counts');
xlim([1 8]);
ylim([0 135]);

subplot(1,3,2);
bar(AxzC,HistContent,'b');
xlabel('cluster DNA content(bp)')
ylabel('counts');
xlim([-1000 4000]);

subplot(1,3,3);
plot(Contents, RadGyr3,'ko', 'MarkerSize',5, 'MarkerFaceColor','y');
xlabel('cluster DNA content(bp)')
ylabel('equivalent diameter, nm');
ylim([0 1000]);







