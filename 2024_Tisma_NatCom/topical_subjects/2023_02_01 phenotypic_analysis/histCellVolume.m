%% Written by Géraldine Laloux, UCLouvain, Feb 2020
% PMID: 34256020
% Edited and used in "Direct observation of a crescent-shape chromosome in Bacillus subtilis" By Tišma et al. 2023
%% Microscope resolution
px2um = 0.065; 

for frame=1:length(cellList.meshData)
    for cell=1:length(cellList.meshData{frame})
        if isfield (cellList.meshData{frame}{cell},'length')
            cellListExtra.meshData{frame}{cell} = cellList.meshData{frame}{cell};
        end
        if ~isfield (cellList.meshData{frame}{cell},'length')
        cellListExtra.meshData{frame}{cell} = GetExtraData(cellList.meshData{frame}{cell});
        end
    end
end
%% gets relevant cell area information for each cell on each frame and creates a vector containing all area values in pixels.
n = 0;
m = 0;
cellL = cellListExtra.meshData;
cellArea = [];
cellVol = [];
for frame=1:length(cellL)
    for cell=1:length(cellL{frame})
        if ~isempty(cellL{frame}{cell}) && isfield(cellL{frame}{cell},'area')
            n = n + 1;
            a = (cellL{frame}{cell}.area)*px2um^2;
            cellArea = [cellArea a];
        end
        if ~isempty(cellL{frame}{cell}) && isfield(cellL{frame}{cell},'volume')
            m = m + 1;
            b = (cellL{frame}{cell}.volume)*px2um^3;
            cellVol = [cellVol b];
        end
    end
end
cellArea = cellArea(~isnan(cellArea));
cellVol = cellVol(~isnan(cellVol));
%% Statistical parameters
CellNumber=n;          % Calculates cell number value. Used for summary table export.
Mean=mean(cellArea); % Calculates mean. Used for summary table export.
Stdev=std(cellArea); % Calculates standard. Used for summary table export.
CV = std(cellArea)/mean(cellArea); % Calculates CV (coefficient of variation). Used for summary table export.
% Volume =  % Calculates CV (coefficient of variation). Used for summary table export.
disp(['cell number:' num2str(n)]);
disp(['meanArea: ' num2str(mean(cellArea))]);
disp(['stdv : ' num2str(std(cellArea))]);
disp(['CV: ' num2str(CV)]);
disp(['mean Volume: ' num2str(mean(cellVol))  ' ± ' num2str(std(cellVol))])
%% creation of histogram data
c = 0:0.5:10; % this gives the range and step size of the histogram (ex: 0:1:15 = from 0 to 15 with steps of 1). Modify it if needed.
h = hist(cellVol,c);
hPERC = 100*h/(sum(h));
%% creation of histogram data
h1 = hist(cellVol1,c);
hPERC1 = 100*h1/(sum(h1));
%% Figure
figure;
hold on;
bar(c,hPERC,'k');
bar(c,hPERC1,'b');
plot(c,hPERC,'k');
plot(c,hPERC1,'b');
ylabel('Fraction of cells (%)','FontSize',20)
xlabel('Cell volume (µm3)','FontSize',20)
axis square
% xlim ([0 1]);  %limits of the X axis; change if needed.
% ylim ([0 50]); % limits of the Y axis; change if needed.
box off
hold off
%% Saving as pdf and fig
saveas(gcf,'name.pdf');
saveas(gcf,'name.fig');