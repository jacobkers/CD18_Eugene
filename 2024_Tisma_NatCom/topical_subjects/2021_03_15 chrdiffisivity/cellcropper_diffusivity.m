clear all;
dirsep='/'; % change it for windows
dircrop0='/Users/fabaiwu/Documents/work/drafts/ChrCompactionSegregation/Figure2/material/diffusivity/20161017shorterm';
dirsave=[dircrop0 dirsep 'allcropped'];
% CS=20; % cropping edge
cd(dircrop0)
dirscan=dir('2016*');
namel=numel(dirscan(1).name);
for scans=1:10
    scanname=dirscan(scans).name;
    dircrop=[dircrop0 dirsep scanname];
    frs=61;
    chs=4; % phase + 3 color channels
    % chbf=1;
    Nstr=1:4; % how many xy positions?
    Npos=1:1; % these two ranges are defined such that you can process part of the data in this folder

    cd(dircrop);
    dirxy=dir('*xy*');
    for i=Nstr
        dirlab1=[dircrop dirsep dirxy(i).name dirsep 'coordinates'];

        cd(dirlab1);
        poses=dir('xy*.mat');
        for j=Npos
            load(poses(j).name);
            for c=1:chs
                cd(dirfiles);
                dircs=dir(cat(2,'*c',int2str(c),'*.tif'));
                for f=1:length(dircs)
                    cd(dirfiles);
                    I=imread(dircs(f).name);
                    for zz=1:size(cropcoord,1)
                        Icrop=I(cropcoord(zz,2):cropcoord(zz,3),cropcoord(zz,4):cropcoord(zz,5));
                        namebase=cat(2,'scan', scanname(namel-2:namel),'_', dirxy(i).name,'_cellnr',num2str(cropcoord(zz,1),'%03i'),'_c',int2str(c),'.tif');
                        cd(dirsave);
                        imwrite(Icrop,namebase,'tiff','Compression','none','WriteMode','Append');
                    end
                end
            end
        end
    end
end
