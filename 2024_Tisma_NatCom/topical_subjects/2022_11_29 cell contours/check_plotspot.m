function occupied=check_plotspot(plotx,ploty,usedplotcenters_xy, radius)
occupied=0; 
if ~isempty(usedplotcenters_xy)
    dist=((usedplotcenters_xy(:,1)-plotx).^2+(usedplotcenters_xy(:,2)-ploty).^2).^0.5;
    sel=find(dist<radius);
    if ~isempty(sel) 
        occupied=1; 
    end
end