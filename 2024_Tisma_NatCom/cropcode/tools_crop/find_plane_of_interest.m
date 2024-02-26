function [plane_of_interest,focal_plane,fp_idx]=find_plane_of_interest(pth, z_name_template,whatplane,sho )
%load full stack to find focal plane 
curpth=pwd;
stck=[]; %raw full range stack from any channel
%set stack planes by hand
cd(pth);
if sho, figure; end
cd(pth);
plane_names=dir(z_name_template);
[Np,~]=size(plane_names);
for plni = 1:Np 
    z_name=plane_names(plni).name;
    newplane=double(imread(z_name));
     stck = cat(3,stck,newplane);
end
%get focal plane  (assuming largely planar features)
stcurve=squeeze(std(std(stck)));
[~,fp_idx]=max(stcurve);
focal_plane=double(stck(:,:,fp_idx));
switch whatplane 
    case 'focalplane'
        plane_of_interest=focal_plane;
    case 'max_project'
        plane_of_interest=squeeze(max(stck,[],3));
end
  
cd(curpth) ; 
if sho
    subplot(1,2,1);     pcolor(focal_plane); shading flat;
    subplot(1,2,2);     pcolor(plane_of_interest); shading flat;
end

    
    