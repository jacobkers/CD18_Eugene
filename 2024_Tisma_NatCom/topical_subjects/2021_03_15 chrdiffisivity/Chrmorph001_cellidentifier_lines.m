%% Note 15/03/2016: this code takes deconvoluted images, crop off single cells, save the labeled data and cell size data
%% This one you can stop and then pick up the analyzed data and keep going again.
%% Analyses on 
% cd('/Users/fabaiwu/Applications/dip')
% run('startup.m');
clear all;
close all;
cmpstr=computer;%code for determining which computer this is running under
if cmpstr(1:3)=='PCW'
    dirsep='\';%for pc
elseif cmpstr(1:3)=='GLN'
    dirsep='/';%for linux
elseif cmpstr(1:3)=='MAC'
    dirsep='/';%for mac (cheack dir seperator)
end
%%%%%% Go to data directory
%dirraw = 'I:\Xuan\20160122\timelapse';
dirraw='O:\Fabai\10minshorterm\20160921bn2179007\xy1\DriftCorrected';
% dirdec='L:\BN\CD\Shared\XuanZheng\20160122\deconvolutedcelldata';
% in both folders the cells are separated into folders t01-t05
% This program does not match cells from t1 to t2, but will record cell
% positions which will allow matching afterwards
dirsave='O:\Fabai\10minshorterm\20160921bn2179007\xy1\crop1';
%dirwrite='L:\BN\CD\Shared\XuanZheng\20160122\Fabaimatlabtest';
%% define several data files
cellindex=[];
f0=1;
p0=1;
label0 = 0;
%% define a few parameters
pxval=15;
CS=16; % cropping space in all directions
MagDisplay=4;
cd(dirraw);
fds=dir('t0*');
for f=1:length(fds); % go through folders
    cd([dirraw dirsep fds(f).name]);
    poses=dir('*xy*');
    fname1=poses(1).name;
    indxy=find(fname1=='y'); % so the number of position is cat(2,'xy',fname1(indxy+1:indxy+2))
    if f==f0;
        posini=p0;
    else
        posini=1;
    end
    for p=posini:length(poses); % go through imaging positions
        cd([dirraw dirsep fds(f).name]);
        pname=poses(p).name;
        xystr=pname(indxy-1:indxy+2);
        xynum=str2double(xystr);
        cd(pname);
        pos4=dir('*c4*');
        numfr=length(pos4);
        %% take the brightfield images and subtract background to produce a inverted image (the cell is now bright)
        % Finding focal plane
        br = zeros(1,3);
        for frame = 14:16
            I1 = imread(pos4(frame).name);
            br(1,frame-13) = max(I1(:));
        end
        focalplane = 13 + find(br(1,:) == max(br(1,:)));
        phase1=double(imread(pos4(focalplane).name)); 
        bg = gaussf(phase1,20,'best');
        Iph = bg - phase1;
        [imw,imh] = size(Iph);
        %% thresholding to separate cells
        thres=0;
        Iph1=Iph>thres;
        lab0=label(Iph1,2,500,10000); % set threshold for the cell size
        lab=int16(lab0);
        maxlab=max(lab(:));
        %% Get rid of cells too close to the edge, simply annoying to crop
        labvals=unique(nonzeros(lab));
        tpbt=lab([1:CS imh-CS:imh],:);lfrt=lab(:,[1:CS imw-CS:imw]);
        edgevals=unique(cat(1,tpbt(:),lfrt(:))); % find labels that are close to the edge;
        Getrid=ismember(labvals,edgevals);
        labstay=find(Getrid==0);
        if numel(labstay)==0; continue; % if all cells are excluded then go to the next loop
        else
            labdel=labvals(Getrid==1);
            for del=1:numel(labdel);
                lab(lab==labdel(del))=0;
            end
            %% now point out the double cells using the mouse (left: select, right: stop selection)
             but=1;
             dublab=[];
             dipshow(lab,'Labels');
             while but==1;
                 [xi,yi,but]=ginput(1);
                 xi=round(xi);yi=round(yi);
                 if but==1
                    dub=lab(yi,xi); % find the label number of the dublet cell
                    if dub==0;
                        dub=mode(double(nonzeros(lab(yi-3:yi+3,xi-3:xi+3)))); % in case the clicking is off
                    end
                    dublab=cat(1,dublab,dub);
                 end
             end
             close(gcf);
             %% now take the dublets out and separate them one by one
             labnew=maxlab;
             dd=1;
             while dd<=numel(dublab);
                 [ys,xs]=find(lab==dublab(dd));
                 cc=[min(ys)-CS max(ys)+CS min(xs)-CS max(xs)+CS]; % cropping coordinates, y=cc(1)+y0-1; x=cc(3)+x0-1;
                 labcrop0=lab(cc(1):cc(2),cc(3):cc(4)); % original, keep
                 labcrop=dip_image(labcrop0);
                 numlab2=1;
                 Ilabcrop=(labcrop==dublab(dd));
                 while numlab2==1;
                     Ilabcrop2=imresize(double(Ilabcrop),MagDisplay,'bilinear');
                     dipshow(dip_image(Ilabcrop2),'Labels');
                     buta=1;
                     linex=[];
                     liney=[];
                     while buta==1;
                          [xi,yi,buta]=ginput(1);
                          if buta==1
                             linex=cat(1,linex,xi); 
                             liney=cat(1,liney,yi);
                          end
                     end
                     close(gcf);
                     linex=linex./MagDisplay;liney=liney./MagDisplay;
                     ydist=max(liney)-min(liney);
                     xdist=max(linex)-min(linex);
                     maxdist=sqrt((ydist)^2+(xdist)^2);
                     if xdist>=ydist;
                         xq=min(linex):(xdist/maxdist)/2:max(linex);
                         vq1=interp1(linex,liney,xq,'spline');
                         vqs=floor([xq;vq1])'; % vq is the [y x] coordinates of all interpolated points
                     else
                         yq=min(liney):(ydist/maxdist)/2:max(liney);
                         vq1=interp1(liney,linex,yq,'spline');
                         vqs=floor([vq1;yq])';
                     end
                     for vs=1:size(vqs,1);
                         Ilabcrop(vqs(vs,1):vqs(vs,1)+1,vqs(vs,2):vqs(vs,2)+1)=0;
                     end
                     lab2=label(Ilabcrop,2,200,10000);
                     numlab2=unique(nonzeros(double(lab2)));
                     Ilabcrop=lab2;
                 end
                 dipshow(lab2,'Labels');
                 %% do a quality control over the new labeling, if not, do it again
                 quest_dualsep = questdlg('Is the data OK?','Quality Control','Yes','Do again','Move on','Yes');
                    switch quest_dualsep
                        case 'Yes'
                            close(gcf);
                            Ilabcrop=int16(Ilabcrop);
                            labnew=labnew+1;
                            Ilabcrop(Ilabcrop==1)=dublab(dd);
                            Ilabcrop(Ilabcrop==2)=labnew;
                            %% Now put the new label back in the original labeled picture
                            labcrop0(labcrop0==dublab(dd))=Ilabcrop(labcrop0==dublab(dd));
                            lab(cc(1):cc(2),cc(3):cc(4))=labcrop0;
                            data_qual = 0;
                            dd=dd+1;
                        case 'Do again'
                            close(gcf);
                        case 'Move on'
                            close(gcf);
                            lab(cc(1):cc(2),cc(3):cc(4))=0;
                            dd=dd+1;
                    end
             end
            %% Here we need to add a step to draw boundaries between objects that are connected
            lab=lab+label0;
            lab(lab==label0)=0;
            label0=label0+labnew;
            data=measure(lab,[],{'size','feret','center','minimum','maximum'}); % minimum/maximum are shown in [x;y]
            % ferets data: maxFeret, minFeret, perpenFeret, maxFangle, minFangle (note, angles are -2pi to 2pi)
            labvals=double(unique(nonzeros(lab)));
            %% later should save all the labeled cells in one folder & their info in one matrix
            celldat1=cat(2,labvals,double(data.size)',double(data.feret)',(data.center)',(data.minimum)',(data.maximum)');
            cellindex=cat(1,cellindex,celldat1);    
            %% Now store the images and the data
%            cd(dirsave);
            fdnum=f;
            posnum=p;
            % save(cat(2,'fd',fds(f).name,xystr,'.mat'),'lab','numfr','xynum','fdnum','celldat1','posnum','label0');
            close(gcf);
        end
        % save('celldata.mat','cellindex','fdnum','posnum');
    end
end

