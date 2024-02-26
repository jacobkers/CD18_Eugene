function pathinfo = X000_setpath4snapshots(expno,usr);
%%% THIS IS FOR SNAPSHOTS DATA WITHOUT DIFFERENT TIME POINTS
%%% FIRST OF ALL: TELL THE CODE WHERE TO GO
%%% 20170313
% This code is used for setting the pathinfo for the dataset
% Reminder 1: Avoid using 't','y','c','z' in the names of files except for
% indicating time points or region of interest
% Reminder 2: This code assumes all the raw data(different
% timepoints,different regions,different channels,different z) are all
% saved in pathinfo.dirraw
% Reminder 3: naming order: t,xy,c,z


switch usr
    case 'Sandro',      pathinfo=X_User_config_Sandro(expno);
    case 'Jacob' ,      pathinfo=X_User_config_Jacob(expno);
    case 'Sandro_2,',   pathinfo=X_User_config_Sandro_2(expno);       
end

addpath(genpath(pathinfo.dircode));
addpath(pathinfo.dirdipimage);
cd ..; addpath(genpath([pwd,'\common_tools\'])); 
cd(pathinfo.dircode);

% if not specified before:
if ~isfield(pathinfo, 'dirraw') 
    pathinfo.dirraw =            pathinfo.experimentpath; % Raw images exported from the microscope
end

% inferred settings:
addpath(genpath(pathinfo.dircode));
pathinfo.dirsep = '\'; % Windows
pathinfo.dirdeconvolution = [pathinfo.experimentpath, 'deconvolve']; % Deconvolved images in fluorescence channels
pathinfo.dircellcoordinate =[pathinfo.experimentpath, 'X020_cellcoordinate'] ; % The coordinates of cells obtained after cell identification
pathinfo.dircrop =          [pathinfo.experimentpath, 'X030_cellcrop']; % The cropped cells
pathinfo.dirgrid =          [pathinfo.experimentpath, 'X040_cellgrid']; % The grids generated
pathinfo.dirmeasurement =   [pathinfo.experimentpath, 'X040_measurement']; % The measurement file
pathinfo.dirdensityanalysis=[pathinfo.experimentpath, 'X050_densitydata']; % Selected cells for ring-chromosome DNA density analysis
pathinfo.xydigit = '%02i'; % how many digits there are for the 'xy' naming
pathinfo.channel = 1:pathinfo.numberofcolours; % Channel order: Phase,chromosome,ori,ter, includes phase
pathinfo.zdigitdeconvolution=repmat(['%03i'],pathinfo.numberofcolours,1); % how many digits there are for the 'z' naming for each channel
%%% DECONVOLUCTION  CHROMOSOME;ORI;TER
pathinfo.timelapse = 0; % 1: Time lapse(including different rois taken at different time points) 0: Snapshots
pathinfo.date = '20171012'; % On which day the analysis was conducted
if ~isfield(pathinfo, 'txycz_template')
        pathinfo.txycz_template=[pathinfo.namebases{1} pathinfo.txycz_base '.tif'];
end



