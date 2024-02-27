%Script for shell shape click work:
% crescent:   less  than ~50% of cell area, one inflexion
% pancake:    larger than 50% of cell area, roughly circular circumference
% multilobe:  multiple inflexions
% compact:    less than 50%, featureless spot
%Note: keep as script to avoid time-consuming re-loading
%Jacob Kerssemakers, 2023/24

actions.re_load=0; %for re-loading oufti-data
actions.re_run=0;  %for re-doing analysis (otherwise, just plotting existing)

%for user-click classification
actions.user_judge=0;       %
actions.re_judge_append=0;  
%auto-overview:
actions.plot_shapes=1;

    %% Microscope resolution (um to pixels conversion) - modify this value when required
    %px2um = 0.07;  Dumbledore 1.5X pixel size
    %px2um = 0.11;  Dumbledore pixel size
    %px2um = 0.04;  Snape 1.5X pixel size
    %px2um = 0.065; Snape pixel size

%% add paths
% paths
codepth=pwd;
addpath(codepth);
cd ..; cd ..;
addpath(genpath([pwd,'\common_tools\'])); 
cd(codepth);

thistopic='C:\Users\jkerssemakers\Dropbox\CD_recent\BN_CD22_Tisma\2022_11_29 cell contours';
close all;

expi=1;
banana_init=AT000_Init_experiment(expi);

if actions.re_load
    load([banana_init.datapath_in banana_init.datafile]);
end

if actions.re_run       
       %init:
        plotinfo.used_xy=[];
        plotinfo.non_overlap=0.15;
        plotinfo.scalar=100;
        counters.gen = 0;
        %check existing data:
        if ~actions.re_judge_append
            shape_data=[];
            counters.click=0;
            counters.last_click_count=0;
        else
            load(banana_init.datapath_out, [banana_init.exp_label, '_shape_data.mat']);
            clicked=find(shape_data(:,3)>0);
            counters.last_click_count=length(shape_data(clicked,1));
            counters.click=last_click_count;
        end  
        switch banana_init.data_type  %walk data, depending on input structure;
            case 'oufti'
            %% 1. This section gets extradata from the meshData provided by Oufti (to obtain calculated length, width, ... among other cell features)
            % and creates a new cellList named cellListExtra
            % Bananalyzer_CellParametersExtract_jk
            % original  ritten by Géraldine Laloux, UCLouvain, Feb 2020
            % expanded by % Jacob Kerssemakers, 2023
            %uses cellList, cellL
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
            cellL = cellListExtra.meshData;
            %% 2. walk the oufti cells to do hand_clicks
            for frame=1:length(cellL)
                for cell=1:length(cellL{frame})
                    if counters.gen<banana_init.total_sampling
                        if ~isempty(cellL{frame}{cell}) && isfield(cellL{frame}{cell}.objects,'area') && isfield(cellL{frame}{cell},'area')
                            %% obtain oufti data:
                            if length(cellL{frame}{cell}.objects.area) == 1
                                counters.gen = counters.gen + 1;    
                                %area:
                                na = cellL{frame}{cell}.objects.area{1};                             
                                %fetch from oufti data:                               
                                chrX=banana_init.px2um*cellList.meshData{1, frame}{1, cell}.objects.outlines{1}(:,1);
                                chrY=banana_init.px2um*cellList.meshData{1, frame}{1, cell}.objects.outlines{1}(:,2);
                                cellX=banana_init.px2um*cellListExtra.meshData{1, frame}{1, cell}.model(:,1);
                               cellY=banana_init.px2um*cellListExtra.meshData{1, frame}{1, cell}.model(:,2);
                                BW=cellList.meshData{1, frame}{1, cell}.objects.masks{1,1};           
                                extra_im=[];
                               [shape_data,counters,plotinfo]=get_symmetry_and_user_info(chrX,chrY,cellX,cellY,BW,extra_im, banana_init,actions,shape_data,counters,plotinfo);
                            end
                        end
                    end
                end
            end
        end
        save([banana_init.datapath_out, banana_init.exp_label, '_shape_data.mat'], 'shape_data'); 
end

%post-processing section:
AT030_Bananalyzer_Post_Process(banana_init);

