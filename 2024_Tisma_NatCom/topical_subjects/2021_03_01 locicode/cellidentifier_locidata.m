% This code takes drift-corrected data of phase contrast images, average them and
% identify cells, then it saves the coordinates and parameters of these
% cells for cropping both phase contrast images and time-lapse images of
% GFP-labeled loci.
clear all;
dirraw = 'E:\Sandro\2017_08_18_2442slmA-rod-dynamics\2442slma-rod-dynamics-001.nd2001\c1'; % Depends on where is raw data is saved
dirsave='D:\My Data\Microscopy\Dynamics\Cropcells';% Depends on where the processed results: mat files will be saved

strain='2762'; % name string to start with.
cellindex=[];
f0=2; % which folder to start now
label0 = 0; % How many cells have been saved  
CS=20; % cropping space in all directions
MagDisplay=4; % To enlarge the cropped region for better visualization and drawing

cd(dirraw);
%%% Go through files found
fds=dir(['*' strain '*']); 
dirsep = '\'; % if Windows computer, change to '\'.
for f=f0:length(fds);
    fstr=num2str(f,'%03i'); % for naming saved files
    dirfiles=([dirraw dirsep fds(f).name dirsep 'DriftCorrected']);
    cd(dirfiles);
    frs = dir('*c1.tif');
    frsloci=dir('*c2.tif');
    Iphall=[];
    for i=1:10;%length(frs);
        Iphall=cat(3,Iphall,double(imread(frs(i).name)));
    end
    Iloci=double(imread(frsloci(1).name));
    phase1=mean(Iphall,3);
    bg = gaussf(phase1,20,'best'); % Using Gaussian to smooth the image to be analyzed
    Iph = bg - phase1;
    [imw,imh] = size(Iph);
    % thresholding to find cells
    thres=0.1*max(Iph(:));
    Iph1=Iph>thres;
    lab0=label(Iph1,2,300,10000); % set threshold for the cell size, smallest size should be 500 otherwise some cells will be excluded
    lab=int16(lab0);
    dipshow(lab,'Labels');
    close(gcf);
    %pause;
    lab1 = lab;
    lab1(lab1>0) = 0;
    maxlab=max(lab(:));
    % Get rid of cells too close to the edge, simply annoying to crop
    labvals=unique(nonzeros(lab));
    tpbt=lab([1:CS imh-CS:imh],:);lfrt=lab(:,[1:CS imw-CS:imw]);
    edgevals=unique(cat(1,tpbt(:),lfrt(:))); % find labels that are close to the edge;
    Getrid=ismember(labvals,edgevals);
    labstay=find(Getrid==0);
    if numel(labstay)==0; continue; % if all cells are excluded then go to the next loop
    else
        labdel=labvals(Getrid==1);
        for del=1:numel(labdel);
            lab(lab==labdel(del))=0; % the exclued cells are labelled as 0
        end
        cellnum = numel(labstay);
        labind = 1;
        label1 = 0;
        while labind <= maxlab
            if isempty(find(lab == labind, 1)) % 0 means not empty
                labind = labind +1;
            else
                [label1,lab1,labind]=cellprescreening_loci(lab,lab1,labind,phase1,CS,label1,MagDisplay,Iloci);
            end
        end
        label0 = label0 + label1;
        % dipshow(lab1,'Labels');
        % pause;
        % close(gcf);
        data=measure(lab1,[],{'size','feret','center','minimum','maximum'}); % minimum/maximum are shown in [x;y]
        % ferets data: maxFeret, minFeret, perpenFeret, maxFangle, minFangle (note, angles are -2pi to 2pi)
        labvals=double(unique(nonzeros(lab1)));
        %% later should save all the labeled cells in one folder & their info in one matrix
        celldat1=cat(2,labvals,double(data.size)',double(data.feret)',(data.center)',(data.minimum)',(data.maximum)');
        cellindex=cat(1,cellindex,celldat1);    
        %% Now store the images and the data
        dircrop=[dirsave dirsep 'BN' strain];
        cd(dircrop);
        save(cat(2,'BN',strain,'_',fstr,'labels.mat'),'lab1','celldat1','label0','dirfiles','frs','frsloci','phase1','strain');
        close(gcf);
    end
    save('celldat1.mat','cellindex');
end



