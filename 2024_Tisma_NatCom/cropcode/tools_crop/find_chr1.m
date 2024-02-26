function [Chrdat,Chrbi,Chrsum,Chrcent] = find_chr1(Ilab);
% note now Ilab is a double 
% Ilab=Chrlab1;
% labthres=1;
% This function assumes that all the labeled objects in the matrix belongs
% one single chromosome, so it will connect different parts and finally
% output an artificially connected chromosome for size measurement.
% However, note that the geometric center and the mass should be determined
% by the original chromosome without artificial connections.
Ibi=Ilab>=1;
%now find the leftmost/rightmost points and connect them, and find
%top/bottom points and connect them
[ys,xs]=find(Ibi==1);
Chrsum=numel(ys);Chrcent=[mean(xs) mean(ys)];
% [ytop,indtop]=max(ys);xtop=xs(indtop);
% [ybot,indbot]=min(ys);xbot=xs(indbot);
% [xleft,indleft]=min(xs);yleft=ys(indleft);
% [xright,indright]=max(xs);yright=ys(indright);
% coordsA=[ytop xtop; yright xright; ybot xbot; yleft xleft];
Ilab2=label(Ibi,2,12,10000);
numlab=max(Ilab2(:));
Ibi=double(Ilab2)>0;
% Labmsr=measure(Ilab2,[],{'center'});
% Labcents=double(Labmsr.center);
% coordsA=Labcents';
coordsA=[];
for j=1:numlab;
    [y1s,x1s]=find(double(Ilab2)==j);
    numpos=round(numel(y1s)/2);
    cent1=[y1s(numpos) x1s(numpos)];
    coordsA=cat(1,coordsA,cent1);
end
if numlab>1;
for i=1:numlab-1;
    ydist=coordsA(i+1,1)-coordsA(i,1);
    xdist=coordsA(i+1,2)-coordsA(i,2);
    maxdist=sqrt((ydist)^2+(xdist)^2);
                     if abs(xdist)>=abs(ydist);
                         xq=coordsA(i,2):(xdist/maxdist)/2:coordsA(i+1,2);
                         vq1=interp1(coordsA([i i+1],2),coordsA([i i+1],1),xq,'spline');
                         vqs=floor([xq;vq1])'; % vq is the [y x] coordinates of all interpolated points
                     else
                         yq=coordsA(i,1):(ydist/maxdist)/2:coordsA(i+1,1);
                         vq1=interp1(coordsA([i i+1],1),coordsA([i i+1],2),yq,'spline');
                         vqs=floor([vq1;yq])';
                     end
                     for vs=1:size(vqs,1);
                         Ibi(vqs(vs,2):vqs(vs,2)+1,vqs(vs,1):vqs(vs,1)+1)=1;
                     end
end
end
Chrbi=Ibi;
Chrdat=measure(Chrbi,[],{'size','feret','center'});
%handle empty measurement:
if ~isfield(Chrdat,'size')
    Chrdat=struct('size', NaN);
    if ~isfield(Chrdat,'size'),Chrdat.size=NaN; end 
    if ~isfield(Chrdat,'feret'),Chrdat.feret=[NaN ; NaN ; NaN ; NaN; NaN]; end 
    if ~isfield(Chrdat,'center'),Chrdat.center=[NaN NaN]; end 
end
% Chrdat.size=Chrsum; Chrdat.center=Chrcent;
% dipshow(Chrbi,'Labels');
% Chrmsr=cat(2,double(Chrsum),double(Chrdat.feret)',double(Chrcent));


