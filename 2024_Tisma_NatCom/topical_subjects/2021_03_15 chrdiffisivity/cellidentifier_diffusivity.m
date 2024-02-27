% This code takes drift-corrected data of phase contrast images, average them and
% identify cells, then it saves the coordinates and parameters of these
% cells for cropping both phase contrast images and time-lapse images of
% GFP-labeled loci.
clear all;
dirsep = '/'; % if Windows computer, change to '\'.
dirraw0 = '/Users/fabaiwu/Documents/work/drafts/ChrCompactionSegregation/Figure2/material/diffusivity/20161017shorterm'; % Depends on where is raw data is saved
% dirsave='D:\fabaiwu\Desktop\chrloci\croppedcells';% Depends on where the processed results: mat files will be saved
cd(dirraw0)
dirpos=dir('2016*');
for p=1:length(dirpos)
    p
    dirraw=[dirraw0 '/' dirpos(p).name];
    posb='xy'; % name string to start with.
    cellindex=[];
    f0=1; % which folder to start now
    label0 = 0; % How many cells have been saved  
    labpick=2; % pick the last frame to process because cells grow
    CS=20; % cropping space in all directions
    MagDisplay=4; % To enlarge the cropped region for better visualization and drawing

    cd(dirraw);
    %%% Go through files found
    fds=dir([posb '*']); 
    for f=f0:length(fds)
        f
        cropcoord=[];
        fstr=num2str(f,'%03i'); % for naming saved files
        dirxy=[dirraw dirsep fds(f).name];
        cd(dirxy);
        mkdir('coordinates');
        dirfiles=[dirxy dirsep 'DriftCorrected'];
        dirsave=[dirxy dirsep 'coordinates'];
        cd(dirfiles);
        frs = dir('*c1*.tif');
        frsloci=dir('*c3*.tif');
        frshu=dir('*c2*.tif');
        ts=length(frs);
        alllabs=[];Iloci=[];Ihu=[];phall=[];
        for i=[1 ts]
            % read files
            phase1=double(imread(frs(i).name));
            phall=cat(3,phall,phase1);
            Iloci=cat(3,Iloci,double(imread(frsloci(i).name)));
            Ihu=cat(3, Ihu, double(imread(frshu(i).name)));
            ph=gaussf(phase1,2,'best');
            bg = gaussf(phase1,10,'best'); % Using Gaussian to smooth the image to be analyzed
            Iph = bg - ph;
            [imw,imh] = size(Iph);
            % thresholding to find cells
            thres=0.3*max(Iph(:));
            Iph1=Iph>thres;
            lab0=label(Iph1,2,300,10000); % set threshold for the cell size, smallest size should be 500 otherwise some cells will be excluded
            lab=int16(lab0);
            dipshow(lab,'Labels');
            %Quality control over the filtered&labeled images
            quest_dualsep = questdlg('are the isolations OK?','Quality Control','Yes','No','Yes');
                    switch quest_dualsep
                      case 'Yes'
                            close(gcf);
                            alllabs=cat(3,alllabs,lab);
                        case 'No'
                            close(gcf);
                            continue
                    end
        end
        lab = alllabs(:,:,labpick);
        % lab1(lab1>0) = 0;
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
            for del=1:numel(labdel)
                lab(lab==labdel(del))=0; % the exclued cells are labelled as 0
            end
            cellnum = numel(labstay);
            labind = 1;
            label1 = 0;
            while labind <= maxlab
                if isempty(find(lab == labind, 1)) % 0 means not empty
                    labind = labind +1;
                else
                    [r,c] = find(lab == labind);
                    bot=min(r(:))-CS; top=max(r(:))+CS; left=min(c(:))-CS; right=max(c(:))+CS;
                    yesorno=cellprescreening_diffusivity(phall,alllabs,Iloci,Ihu,bot,top,left,right);

                    if yesorno==1
                        cropcoord=cat(1,cropcoord,[labind bot top left right]);
                    else
                        lab(lab==labind)=0;
                    end
                    labind=labind+1;
                end
            end
            label0 = label0 + label1;

            data=measure(lab,[],{'size','feret','center','minimum','maximum'}); % minimum/maximum are shown in [x;y]
            % ferets data: maxFeret, minFeret, perpenFeret, maxFangle, minFangle (note, angles are -2pi to 2pi)
            labvals=double(unique(nonzeros(lab)));
            %% later should save all the labeled cells in one folder & their info in one matrix
            celldat1=cat(2,labvals,double(data.size)',double(data.feret)',(data.center)',(data.minimum)',(data.maximum)');
            cellindex=cat(1,cellindex,celldat1);
            %% Now store the images and the data
            % dircrop=[dirsave dirsep strain];
            cd(dirsave);
            save(cat(2,posb,'_',fstr,'labels.mat'),'lab','celldat1','label0','dirfiles','cropcoord');
            close(gcf);
            % dipshow(lab,'Labels');
        end
        save('celldat1.mat','cellindex');
    end
end