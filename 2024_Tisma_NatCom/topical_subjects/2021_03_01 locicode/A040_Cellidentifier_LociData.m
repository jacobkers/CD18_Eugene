function A040_Cellidentifier_LociData(initval)
%JWJK_A:-------------------------------------------------------------------
%
%Cellidentifier
%
%Summary: code takes drift-corrected data of phase contrast images, average them and
% identify cells, then it saves the coordinates and parameters of these
% cells for cropping both phase contrast images and time-lapse images of
% GFP-labeled loci.
%
%Approach: first, a labeled cell area image is pre-screened by user; 
%areas are accepted, trashed or split when necessary. Then the cleaned
%lalbel image is used for cell measurements. Cel area data is saved for
%later cropping (A50)
%necessary
%
%Input: drift-corrected cell images of supposed movie frame, containing color channels 
%
%Output: saved .mat files for later use. BN_labels.mat contains:
    % label1: updated counter
    %lab1: re-labeled label picture, including separated cells
    %celldat1: some measures of cell size
    %label0: number of cells saved
    % imindir: file path infiles
    % frs: list of brightfield input files
    % frsloci_cfp: list of cfp inputfiles
    % phase1;
%
%References: written by F.WU, edited by SJ/JWJK
%
%:JWJK_A-------------------------------------------------------------------
%dip init
if 1
addpath('C:\Program Files\DIPimage 2.9\common\dipimage');
dip_initialise;
dipsetpref('ImageFilePath', 'C:\Program Files\DIPimage 2.8\images');
end
curpath=pwd;
imindir =[initval.imagepath 'A20_DriftCorrected\'];



cellindex=[];
f0=1; % which folder to start now
label0 = 0; % How many cells have been saved  
CS=initval.cropping_edge; % cropping space in all directions
MagDisplay=4; % To enlarge the cropped region for better visualization and drawing

channo=length(initval.SaveLabels);


dirsep = '\'; % if Windows computer, change to '\'.

%% build averaged phase or brightfield image; perform tresholding to find separate objects
cd(imindir); 

frs = dir(strcat('*', char(initval.SaveLabels{1}),'.tif'));
%list of images per channel (can be empty)
frsloci_color1=dir(strcat('*', char(initval.SaveLabels{2}),'.tif')); 
frsloci_color2=dir(strcat('*', char(initval.SaveLabels{3}),'.tif')); 
frsloci_color3=dir(strcat('*', char(initval.SaveLabels{4}),'.tif')); 
Nframes=min([length(frs), initval.shortset]);


Iphall=[];
for ii=1:Nframes;    %phase or BF images 
    Iphall=cat(3,Iphall,double(imread(frs(ii).name)));
end
cd(initval.codepth);
phase1=mean(Iphall,3);  %average phase or BF image
bg = gaussf(phase1,20,'best'); % Using Gaussian to smooth the image to be analyzed
Iph = bg - phase1;
[imw,imh] = size(Iph);
% thresholding to find cells
thres=0.1*max(Iph(:));
Iph1=Iph>thres;
lab0=label(Iph1,2,300,10000); % set threshold for the cell size, 
%smallest size should be 500 otherwise some cells will be excluded

%lab0 is phase or BF image; links to bflabel; savelalbel(1)
%Iloci is savelabel(2)  (fluorescence);


lab=int16(lab0); 
dipshow(lab,'Labels');
close(gcf);
%% use
cd(imindir); 



Iloci=double(imread(frsloci_color1(1).name));
cd(initval.codepth);
%pause;
lab1 = lab;
lab1(lab1>0) = 0;  %empty picture?
maxlab=max(lab(:));
% Get rid of cells too close to the edge, simply annoying to crop
labvals=unique(nonzeros(lab));  %levels used for labeled cells
tpbt=lab([1:CS imh-CS:imh],:);lfrt=lab(:,[1:CS imw-CS:imw]);
edgevals=unique(cat(1,tpbt(:),lfrt(:))); % find labels that are close to the edge;
Getrid=ismember(labvals,edgevals);
labstay=find(Getrid==0);
labdel=labvals(Getrid==1);
for del=1:numel(labdel);
    lab(lab==labdel(del))=0; % the exclued cells are labelled as 0
end
cellnum = numel(labstay);
labind = 1;        
%label index
label1 = 0;
while labind <= maxlab  %work all labels for this image
    disp(strcat('Cell',num2str(labind),'of',num2str(maxlab)));
    if isempty(find(lab == labind, 1)) % 0 means not empty
        labind = labind +1;
    else
        [label1,lab1,labind]=cellprescreening_loci(lab,lab1,labind,phase1,CS,label1,MagDisplay,Iloci,initval);
    end
end
labind = maxlab;
label0 = label0 + label1;

%show all the detected labels
dipshow(lab1,'Labels');
pause(3);
close(gcf);

%here, data is measured on the cell (averaged BF or Phase image)
data=measure(lab1,[],{'size','feret','center','minimum','maximum'}); % minimum/maximum are shown in [x;y]
% ferets data: maxFeret, minFeret, perpenFeret, maxFangle, minFangle (note, angles are -2pi to 2pi)

labvals=double(unique(nonzeros(lab1)));
%% later should save all the labeled cells in one folder & their info in one matrix
celldat1=cat(2,labvals,double(data.size)',double(data.feret)',(data.center)',(data.minimum)',(data.maximum)');
cellindex=cat(1,cellindex,celldat1);    
%% Now store the images and the data
testname='frs';
save(strcat(initval.matfilesinoutdir,'A040_BN','_','labels.mat'),'lab1','celldat1','label0',...
                       'imindir','frs','frsloci_color1','frsloci_color2','frsloci_color3','phase1');
save(strcat(initval.matfilesinoutdir,'A040_celldat1.mat'),'cellindex');
close(gcf);
cd(curpath);



