
%% From the original publication of Kaljevic et al. 2021 (PMID: 34256020)

% Edited by M. Tisma and J. Kerssemakers and used in "Direct observation of a crescent-shape chromosome in Bacillus subtilis" By Tišma et al. 2023

%% add paths
addpath(genpath('File_path\'));
load('spot_detection_file.mat');

% INPUTS
scalefactor = 0.065; % conversion factor from pix to um
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
%% renames cellListExtra.meshData
cellL = cellListExtra.meshData;
%% script
dArray = [];
n = 0; % cell counter
for frame = 1:length(cellL)
    for cell = 1:length(cellL{frame})
           if    ~isempty (cellL{frame}{cell}) && isfield(cellL{frame}{cell},'spots') && isfield(cellL{frame}{cell},'spots2')...
                 && ~isempty(cellL{frame}{cell}.spots.positions) && ~isempty(cellL{frame}{cell}.spots2.positions)... % usual conditions
                 && length(cellL{frame}{cell}.spots.positions) >= 1 && length(cellL{frame}{cell}.spots2.positions) >= 1  % number of spots signal 1 and 2 must be more than 1
            n = n+1;
            X = [];
            Y = [];
            X2 = [];
            Y2 = [];
            min_12_d = [];
            min_12_i = [];
            N1=length(cellL{frame}{cell}.spots.x);
            N2=length(cellL{frame}{cell}.spots2.x);
            for i_1=1:N1
                %% gets euclidian coordinates of the spot of signal 1
                    X = cellL{frame}{cell}.spots.x(i_1);
                    Y = cellL{frame}{cell}.spots.y(i_1); 
                    dd_12=[];
                    for i_2=1:N2
                        %% gets euclidian coordinates of each spot of signal 2
                            X2 = cellL{frame}{cell}.spots2.x(i_2);
                            Y2 = cellL{frame}{cell}.spots2.y(i_2);
                        %% gets the distances from spot 1 to all spots 2. 
                            dd_12(i_2) = sqrt(((X-X2)^2)+((Y-Y2)^2));       
                    end
                    %nearest spot 2 distance &  index
                    [d_12_min, i_12_min]=min(dd_12);
                    min_12_d(i_1) = d_12_min;
                    min_12_i(i_1)=  i_12_min;                     
            end
             %check for double pointers, discard longer distances
             min_12_d_nw = 0*min_12_d;
             min_12_i_nw = 0*min_12_i;
             for i_1=1:N1
                 this_i2=min_12_i(i_1);
                 this_d2=min_12_d(i_1);
                 similar_i2=find(this_i2==min_12_i);
                 dist_to_these=min_12_d(similar_i2);
                 [dmin,ixmin]=min(dist_to_these);
                 %keep only if equal to shortest:
                 if dmin==this_d2
                     min_12_d_nw(i_1)=dmin;
                     min_12_i_nw(i_1)=this_i2;
                 else
                     min_12_d_nw(i_1)=NaN;
                     min_12_i_nw(i_1)=NaN;
                 end                            
             end
             %add results:
             dArray = [dArray scalefactor*min_12_d_nw];
           end
    end
end
%clean it:
dArray=dArray(~isnan(dArray));
%% OUPUTS
%% plots histogram of distances (in um) between the signal 1 spot and the signal 2 spot
c = 0:0.1:2;

h = hist(dArray,c);
hperc = 100*h/sum(h);
figure;
bar(c,hperc, 0.7);
xlim([-0.2,2])

axis square
%% display total number of cells used for calculation; and total number of spots of signal 1 taken into account.
disp(['number of cells =' num2str(n)])