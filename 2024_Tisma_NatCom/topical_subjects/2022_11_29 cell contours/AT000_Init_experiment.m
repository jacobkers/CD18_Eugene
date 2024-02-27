function init=A000_Init_experiment(expi)

init.total_sampling =1400;  %to do all, use number> number of cells
init.user_sampling=1400;
init.px2um = 0.065; 

switch expi
    case 1
        %settings
        init.exp_label='BSG217_oufti';
        init.data_type='oufti';
        init.datapath_in='BSG217\';
        init.datapath_out='BSG217\';
        init.datafile='BSG217_stack_Outlines_Corr_ObjDetection.mat';
        addpath(genpath('M:\tnw\bn\cd\Shared\Tisma\Bacillus subtilis chromosome imaging\scripts\oufti_windows\source'));
        %settings
        init.blowup=10;  %over-sample images for smoother contours
        init.total_sampling =1400;  %to do all, use number> number of cells
        init.user_sampling=1400;
        init.px2um = 0.07;
        %% Microscope resolution (um to pixels conversion) - modify this value when required
        %px2um = 0.07;  Dumbledore 1.5X pixel size
        %px2um = 0.11;  Dumbledore pixel size
        %px2um = 0.04;  Snape 1.5X pixel size
        %px2um = 0.065; Snape pixel size

end
