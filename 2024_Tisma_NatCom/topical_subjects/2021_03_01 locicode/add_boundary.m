function [Ilabcrop1] = add_boundary(Ilabcrop,MagDisplay)
%JWJK_B:-------------------------------------------------------------------
%User-stet boundary between adjacent cells
%
%Summary: this function allows a user to sepearate two connencted cells by
%clicking boundaries between them
%
%Approach
%1) rightclick: right-click first target cell
%2) Left-click a line to cut first cell; rightclick to end
%3) rightclick: right-click second target cell 
%4) Left-click a line to cut 2nd cell; rightclick to end
%
%Input
%
%Output
%
%References
%
%:JWJK_B-------------------------------------------------------------------
Ilabcropforcell1 = Ilabcrop;
Ilabcropforcell2 = Ilabcrop;
% Show the image; enlarged for conovenience
Ilabcroplarge=imresize(double(Ilabcrop),MagDisplay,'bilinear');
dipshow(dip_image(Ilabcroplarge),'Labels');
%% 1)left-click first target cell 
disp('Left-click first target cell ');
[xi,yi,butb]=ginput(1);
xcell1 = round(xi./MagDisplay);
ycell1 = round(yi./MagDisplay);

close(gcf);
dipshow(dip_image(Ilabcroplarge),'Labels');
buta=1;
linex=[];
liney=[];

%% 2) Left-click a line over cells; rightclick to end
disp('Left-click a line over cells; rightclick to end');
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
 Ilabcropforcell1(vqs(vs,1):vqs(vs,1)+1,vqs(vs,2):vqs(vs,2)+1)=0;
end
lab1=label(Ilabcropforcell1,2,200,10000);lab1 = uint16(lab1);
cell1 = lab1(ycell1,xcell1);


%%3) left-click second target cell 
disp('Left-click second target cell ');
dipshow(dip_image(Ilabcroplarge),'Labels');
[xi,yi,butb]=ginput(1);
xcell2 = round(xi./MagDisplay);
ycell2 = round(yi./MagDisplay);

close(gcf);
dipshow(dip_image(Ilabcroplarge),'Labels');
buta=1;
linex=[];
liney=[];


%%4) Left-click a line to cut 2nd cell; rightclick to end
disp('');
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
 Ilabcropforcell2(vqs(vs,1):vqs(vs,1)+1,vqs(vs,2):vqs(vs,2)+1)=0;
end
close(gcf);
lab2=label(Ilabcropforcell2,2,200,10000);lab2 = uint16(lab2);
cell2 = lab2(ycell2,xcell2);
[sx,sy] = size(Ilabcrop);
Ilabcrop1 = zeros(sy,sx);
Ilabcrop1(find(lab1 == cell1)) = 1;
Ilabcrop1(find(lab2 == cell2)) = 2;
dipshow(Ilabcrop1,'Labels');
end