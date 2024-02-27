

close all;
stack_ori=double(imread('marker1_ori.tif'));
im_bf=double(imread('marker1_bf.tif'));
roidata=readtable('Markers.csv');

%corner point:
pos_xx=roidata.X;  
pos_yy=roidata.Y;  %
N_rois=length(pos_xx);
xi=pos_xx(1);
yi=pos_yy(1);




opts.SelectedVariableNames = [1:5]; 
opts.DataRange = '2:11';
xy_data = readmatrix('Spot Tracker Results_marker1.csv')
xx=xy_data(:,2)-xi;
yy=xy_data(:,3)-yi;
subplot(1,2,1); 
    pcolor(im_bf); shading flat; colormap bone; axis equal; axis tight; hold on
    title('BF');
    plot(xx(1),yy(1), 'ro-', 'MarkerFaceColor', 'r');
    plot(xx,yy, 'ro-');
subplot(1,2,2);  
    pcolor(stack_ori(:,:,1)); shading flat; colormap bone; axis equal; axis tight; hold on
    title('Ori');
    plot(xx(1),yy(1), 'bo-', 'MarkerFaceColor', 'b');
    plot(xx,yy, 'bo-');
dum=1