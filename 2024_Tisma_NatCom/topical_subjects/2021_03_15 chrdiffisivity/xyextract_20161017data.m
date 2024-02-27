% This code takes data of each cell and calculate the coordinate of the
% position of the object along the long and short axis based on the angle
% difference between the cell-center-object-center line and the axes lines.
% The data shall be stored in a new list of files, one file per cell
dircrop='/Users/fabaiwu/Documents/work/drafts/ChrCompactionSegregation/Figure2/material/diffusivity/20161017shorterm/allcropped';

cd(dircrop);
datafiles=dir('*data.mat');
for i=1:length(datafiles);
    load(datafiles(i).name);
    newname=[datafiles(i).name(1:22) 'xyco.mat'];
%    plotname=[datafiles(i).name(1:22) 'plot.tif'];
    ts=size(celldata,1);
    if ts==61;
        %% data for chromosome center
        celll=median(celldata(:,9));cellw=median(celldata(:,10));
        chrx0=chrdata(:,3)-celldata(:,3);
        chry0=chrdata(:,4)-celldata(:,4);
        chrangle0=atan(chry0./chrx0);
        ind=find(chrx0<0);
        chrangle0(ind)=chrangle0(ind)+pi;
        chrangle1=chrangle0-celldata(1,12);% angle difference
        chrdist=sqrt((chrx0).^2+(chry0).^2); % distance from cell center
        %% New coordinates assuming the cells are rotated to align their long axis horizontally
        chrx1=chrdist.*cos(chrangle1); 
        chry1=chrdist.*sin(chrangle1);
        
        %% data for Ori
        orix0=focidata(:,3,1)-celldata(:,3);
        oriy0=focidata(:,4,1)-celldata(:,4);
        oriangle0=atan(oriy0./orix0);
        ind=find(orix0<0);
        oriangle0(ind)=oriangle0(ind)+pi;
        oriangle1=oriangle0-celldata(1,12);
        oridist=sqrt((orix0).^2+(oriy0).^2);
        orix1=oridist.*cos(oriangle1);
        oriy1=oridist.*sin(oriangle1);
        
        %% data for ter
        terx0=focidata(:,3,2)-celldata(:,3);
        tery0=focidata(:,4,2)-celldata(:,4);
        terangle0=atan(tery0./terx0);
        ind=find(terx0<0);
        terangle0(ind)=terangle0(ind)+pi;
        terangle1=terangle0-celldata(1,12);
        terdist=sqrt((terx0).^2+(tery0).^2);
        terx1=terdist.*cos(terangle1);
        tery1=terdist.*sin(terangle1);
        save(newname,'celldata','chrdata','focidata','chrx1','chry1','orix1','oriy1','terx1','tery1','celll','cellw');
    end
end
