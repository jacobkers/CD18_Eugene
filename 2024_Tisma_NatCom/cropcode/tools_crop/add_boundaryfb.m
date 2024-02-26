function [addlab,labB] = add_boundaryfb(Ilab2bit,boundtype,labind,addlab,MagDisplay)

linex=[];liney=[];
buta=1; 
while buta==1
  [xi,yi,buta]=ginput(1);
  if buta==1
     linex=cat(1,linex,xi); 
     liney=cat(1,liney,yi);
  end
end
linex=linex./MagDisplay;liney=liney./MagDisplay;
ydist=max(liney)-min(liney);
xdist=max(linex)-min(linex);
maxdist=sqrt((ydist)^2+(xdist)^2);
if xdist>=ydist
    [~,ia,~]=unique(linex);
    linex=linex(ia);liney=liney(ia); % to avoid repeating axis values
 xq=min(linex):(xdist/maxdist)/2:max(linex);
 vq1=interp1(linex,liney,xq,'spline');
 vqs=floor([xq;vq1])'; % vq is the [y x] coordinates of all interpolated points
else
 [~,ia,~]=unique(liney);
 linex=linex(ia);liney=liney(ia); % to avoid repeating axis values
 yq=min(liney):(ydist/maxdist)/2:max(liney);
 vq1=interp1(liney,linex,yq,'spline');
 vqs=floor([vq1;yq])';
end
for vs=1:size(vqs,1)
 Ilab2bit(vqs(vs,2):vqs(vs,2)+1,vqs(vs,1):vqs(vs,1)+1)=boundtype; % note, switched for non-labled data
end
Ilab2bit=imfill(Ilab2bit,'holes');
lab1=label(Ilab2bit,2,200,10000);lab1 = uint16(lab1);
A=unique(lab1);
A=nonzeros(A);
labB=lab1;
labB(lab1==1)=labind;
if A>1
    for i=2:A
        labB(lab1==i)=addlab+i-2;
        addlab=addlab+1;
    end
end