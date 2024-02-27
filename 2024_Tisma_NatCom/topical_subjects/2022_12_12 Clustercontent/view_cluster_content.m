function view_cluster_content
%Load a cluster report table. Each row entry is a single spot that forms a
%component of a cluster. Note: cell labeles may not be unique.....
pth1='M:\tnw\bn\cd\Shared\Jacob\TESTdata_out\2022_Tisma\BatchanalysisResults\';
switch 3
    case 1  %original MS fig S8B/C
        pth2='Results_Tisma_20221009_4623_c3_SingleBranch_Flipped\';
        fle='Tisma_20221009_4623_A050_Cell_ClusterReport_';    
    case 2
        pth2='Results_Tisma_20220802_4595_SyG_c3_SingleBranch_Flipped\';
        fle='Tisma_20220802_4595_SyG_A050_Cell_ClusterReport_';  
   case 3
        pth2='Results_Tisma_20220803_4595_SyG_c3_SingleBranch_Flipped\';
        fle='Tisma_20220803_4595_SyG_A050_Cell_ClusterReport_';  
end

c2=load([pth1,pth2,fle, 'c2.mat'], 'AllCellsAllClusters');
c3=load([pth1,pth2,fle, 'c3.mat'], 'AllCellsAllClusters');
c4=load([pth1,pth2,fle, 'c4.mat'], 'AllCellsAllClusters');

headers_123=[{'index'};	{'label'};	{'cluster number'}];
headers_per_cluster=[
{'content'};
{'diameter_2Rg'};
{'diameter_ContourBW_Equiv'};
{'diameter_ContourBW_2Rg'};
{'area'};
{'density'};
{'xpos'};
{'ypos'};
{'1Ddistpos'};
{'1DBPpos'};
{'Psf_used'}];

cluster_rank=1;
maxspots_c2=[];
maxspots_c3=[];
maxspots_c4=[]

[N_rows2,~]=size(c2.AllCellsAllClusters);
for ci=1:N_rows2-1
    cell_pointer=ci+1;
    col_pointer=length(headers_123)+(cluster_rank-1)*length(headers_per_cluster)+1;
    maxspots_c2(ci)=100*max(c2.AllCellsAllClusters(cell_pointer,col_pointer));
end

[N_rows3,~]=size(c3.AllCellsAllClusters);
for ci=1:N_rows3-1
    cell_pointer=ci+1;
    col_pointer=length(headers_123)+(cluster_rank-1)*length(headers_per_cluster)+1;
    maxspots_c3(ci)=100*max(c3.AllCellsAllClusters(cell_pointer,col_pointer));
end

[N_rows4,~]=size(c4.AllCellsAllClusters);
for ci=1:N_rows4-1
    cell_pointer=ci+1;
    col_pointer=length(headers_123)+(cluster_rank-1)*length(headers_per_cluster)+1;
    maxspots_c4(ci)=100*max(c4.AllCellsAllClusters(cell_pointer,col_pointer));
end

Ncells=length(maxspots_c2);
binax=[0:5:100];
hist_c2=hist(maxspots_c2, binax)/(N_rows2-1)*100;
hist_c3=hist(maxspots_c3, binax)/(N_rows3-1)*100;
hist_c4=hist(maxspots_c4, binax)/(N_rows4-1)*100;



%% build excel:
%build excel table of barplots
savename=[fle, '_clustercontent.xls'];
data_wrapup1=[binax' hist_c2' hist_c3' hist_c4'];

headers_wrapup1=[{'fluorecence, %'},{'DNA cluster fraction (%)'},{'SMC cluster fraction (%)'}, {'ParB cluster fraction (%)'}];
xlswrite(savename, {'number of cells:'}, 'cluster contents', 'A1');
xlswrite(savename, Ncells, 'cluster contents', 'B1');
xlswrite(savename, headers_wrapup1, 'cluster contents', 'A2');
xlswrite(savename,data_wrapup1, 'cluster contents', 'A3');


%% plots
close all;
subplot(1,3,1); bar(binax, hist_c2);
    xlabel('percentage in main cluster');
    ylabel('percentage');
    legend('c2', 'Location', 'NorthOutside')
subplot(1,3,2); bar(binax, hist_c3, 'r');
    xlabel('percentage in main cluster');
    ylabel('percentage');
    legend('c3', 'Location', 'NorthOutside')
subplot(1,3,3); bar(binax, hist_c4, 'm');
    xlabel('percentage in main cluster');
    ylabel('percentage');
    legend('c4', 'Location', 'NorthOutside')
    
jpgname2=[fle, '_clustercontent.jpg'];
saveas(gcf,jpgname2);

plot2svg([fle, '_clustercontent.svg'], gcf);

