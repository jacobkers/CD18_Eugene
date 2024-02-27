function A005_ClickDriftMarkers(initval)
%JWJK_A:-------------------------------------------------------------------
%Determine sample drift in a series of images via a fixed marker (or blob)
%
%Summary: This code loads from a movie directory, allows the user to click 
%a region with a blob and tracks a vector from it.
%
%Approach: User clicks five static objects. Objects should be more or less
%round for good tracking. Position tracking is done by 2D correlation via
%Fourier domain, i.e. images are correlated with their center-point
%mirrored self. Tracking results from four objects are averaged.
%
%Input: images
%
%Output: a txt file 'drift vector' saved to the same folder as the images
%written by JWJK, 2016/21
%
%:JWJK_A-------------------------------------------------------------------

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

%Jacob Kerssemakers 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; 
curpath=pwd;

%% click it
cd(initval.imagepath);                       %read the proper file names in this directory
FilNames=dir('*c1.tif');
cd(curpath);
firstim=double(imread(strcat(initval.imagepath,FilNames(1).name)));
firstim=firstim(:,:,1);
[wo,ho]=size(firstim); frs=length(FilNames);
P_Color(log10(firstim),wo,ho,'hot'); axis equal;
title('log image-Click five stable markers')


[rclick,cclick]=ginput(5);     %pick five ROIs
%save these
cd(initval.imagepath);                       %read the proper file names in this directory
save('driftmarkers.mat','rclick','cclick');
cd(curpath);
end
close(gcf);
disp('done');

