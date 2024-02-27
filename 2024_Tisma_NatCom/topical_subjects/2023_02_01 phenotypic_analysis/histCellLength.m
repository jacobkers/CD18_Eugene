%% Written by Géraldine Laloux, UCLouvain, Feb 2020
% PMID: 34256020

% Edited and used in "Direct observation of a crescent-shape chromosome in Bacillus subtilis" By Tišma et al. 2023
%% Microscope resolution
px2um = 0.065; 
%% This section gets extradata from the meshData provided by Oufti (to obtain calculated length, width, ... among other cell features)
% and creates a new cellList named cellListExtra
for frame=1:length(cellList.meshData)
    for cell=1:length(cellList.meshData{frame})
        if isfield (cellList.meshData{frame}{cell},'length')
            cellListExtra.meshData{frame}{cell} = cellList.meshData{frame}{cell};
        end
        if ~isfield (cellList.meshData{frame}{cell},'length')
        cellListExtra.meshData{frame}{cell} = getextradata(cellList.meshData{frame}{cell});
        end 
    end
end
%% gets relevant cell length information for each cell on each frame and creates a vector containing all length values in pixels.
n = 0;
cellL = cellListExtra.meshData;
cellLength = [];
for frame=1:length(cellL)
    for cell=1:length(cellL{frame})
        if ~isempty(cellL{frame}{cell}) && isfield(cellL{frame}{cell},'length')
            n = n + 1;
            l = px2um*(cellL{frame}{cell}.length);
            cellLength = [cellLength l];
        end
    end
end
%% Statistical parameters
CellNumber=n;          % Calculates cell number value. Used for summary table export.
Mean=mean(cellLength); % Calculates mean. Used for summary table export.
Stdev=std(cellLength); % Calculates standard. Used for summary table export.
CV = (std(cellLength)/mean(cellLength)); % Calculates CV (coefficient of variation). Used for summary table export.
%To display stat
disp(['cell number:' num2str(n)]);
disp(['mean:' num2str(mean(cellLength))]);
disp(['stdv:' num2str(std(cellLength))]);
disp(['CV: ' num2str(CV)]);
%% creation of histogram data
c = 0:0.5:10; % this gives the range and step size of the histogram (ex: 0:1:15 = from 0 to 15 with steps of 1). Modify it if needed.
h = hist(cellLength,c);
hPERC = 100*h/(sum(h));
%% bar graph settings
figure;
plot(c,hPERC,'k')
ylabel('Fraction of cells (%)','FontSize',20)
xlabel('Cell length (µm)','FontSize',20)
axis square
% xlim ([0 2.5]); %limits of the Y axis; change if needed.
% ylim ([0 50]);  %limits of the X axis; change if needed.
box off