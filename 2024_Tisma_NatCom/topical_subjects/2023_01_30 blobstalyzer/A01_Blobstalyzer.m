function A01_Blobstalyzer
%Simple Clicktool for chromosome classification
%JacobKers2023
N_samples=100;
usr=char(inputdlg('who are you?')); 
if ~isdir(usr), mkdir(usr);end
close all;
%curpth=pwd; for ii=1:2, cd ..; end; 
%addpath(genpath('common_tools'));
%cd(curpth);
%read images
savelabel= ['_', usr,'_',datestr(now, 'mmmm_dd_yyyy_HH_MM')];
stack=[];
stack(:,:,1)=double(imread('1.tif'));
stack(:,:,2)=double(imread('2.tif'));
stack(:,:,3)=double(imread('3.tif'));
stack(:,:,4)=double(imread('4.tif'));

blob_pick_list_ori=[];
for ii=1:4
    [pic1_out,thr]=Find_treshold_MD_V2020(stack(:,:,ii));
    BW=(pic1_out>0);
    bwstruct=bwconncomp(BW,8);    %finds 8-fold connected regions.
    stackdata(ii).blobprops=regionprops(bwstruct,...
        'Centroid', 'Area','MajorAxisLength',...
        'MinorAxisLength','Eccentricity','PixelIdxList','BoundingBox'); 
    labelmat=labelmatrix(bwstruct); %label all the points in regions with the region nr
    N_blobs=double(max(labelmat(:)));
    %[image blob used]
    blob_pick_list_ori=[blob_pick_list_ori; [zeros(N_blobs,1)+ii (1:N_blobs)'] zeros(N_blobs,1)];
end
tic
N_blobs_all=length(blob_pick_list_ori(:,1));
N_blobs_left=N_blobs_all;
blobcounter=0;

blobclass=zeros(N_samples,5);
blob_pick_list_to_use=blob_pick_list_ori;
while blobcounter<N_samples && N_blobs_left>0
    notused=find(blob_pick_list_to_use(:,3)==0);
    N_blobs_left=length(notused);
    if N_blobs_left >0
        blob_pick_list_to_use=blob_pick_list_to_use(notused,:);    
        dice=ceil(N_blobs_left*rand(1)); %pick a random blob
        blobcounter=blobcounter+1;
        imno=blob_pick_list_ori(dice,1);
        blobno=blob_pick_list_to_use(dice,2);
        LTRB=stackdata(imno).blobprops(blobno).BoundingBox;
        [thisblob,~,~]=get_blob_box(stack(:,:,imno),LTRB,10);
        pcolor(thisblob); colormap(jet); shading flat;
        title( ['blob:', num2str(blobcounter),'of:', num2str(N_samples)]);
        [rr,cc]=size(thisblob);
        clr='w';
        shft=3;
        text(shft,shft,'1-donut', 'Color', clr, 'FontSize',15);
        text(shft,rr-shft,'2-crescent', 'Color', clr,'FontSize',15);
        text(cc-2*shft,rr-shft,'3-compact', 'Color', clr,'FontSize',15);
        text(cc-2*shft,shft,'4-other', 'Color', clr,'FontSize',15); 
        pause(0.01);
        tic
        switch 1
            case 1 %#clicking      
                [xn,yn,but]=ginput(1);
                clicktime=toc;
                blobtype=check_LTRB(xn,yn,rr,cc);
           case 2 %#typing      
                blobtype=input('enter type'); 
                if ~ismember(blobtype, [1 2 3 4]), blobtype=4; end    
        end
        clicktime=toc;        
        blobclass(blobcounter,:)=[blobcounter blobtype,imno,blobno, clicktime];
        blob_pick_list_to_use(dice,3)=1;  %used!  
    end 
   
    
    
end
 headers=  [{'index'},{'blobtype'},{'imno'},{'blobno'},{'clicktime'}];
 typecodes=[{'1-donut'},{'2-crescent'},{'3-compact'},{'4-undefined'}];
 save([usr,'\blobstalysis',savelabel,'.mat'],'typecodes','blobclass','headers');
 disp(['clicktime:,', num2str(mean(blobclass(:,5)))]);

function blobtype=check_LTRB(xn,yn,rr,cc)
    corners_x=[1 1 cc cc]; %LB-->clockwise
    corners_y=[1 rr rr 1]; 
    dist=((corners_x-xn).^2+(corners_y-yn).^2).^0.5;
    [mindist,blobtype]=min(dist);

function [snapshot,left,top]=get_blob_box(map,LRTB,pad,offset);
%get a nice area around a pore plus original corner coordinate
%option offset allows for picking an area outside this bounding box,
%towards the center of the image

if nargin<4, offset=0; end

[rr,cc]=size(map);
left=max([floor(LRTB(1))-pad;       1]);
right=min([left+ceil(LRTB(3))+2*pad;  cc]);
top=max([floor(LRTB(2))-pad;        1]);
bot=min([top+ceil(LRTB(4))+2*pad;     rr]);


if ~offset
    snapshot=map(top:bot,left:right);
else
    %check direction:
    rightofmid=sign((left+right/2)-rr/2);  %+1:right of center
    aboveofmid=sign((top+bot/2)-cc/2);  %+1:left of center
    offset_hor=round(-rightofmid*((right-left)+pad));
    offset_ver=round(-aboveofmid*((bot-top)+pad));
    loy=top+offset_ver;
    hiy=bot+offset_ver;
    lox=left+offset_hor;
    hix=right+offset_hor;
    isokay=(lox>=1&hix<=cc&loy>=1&hiy<=rr);
    if isokay
        snapshot=map(loy:hiy,lox:hix);
    else
        %snapshot=0*map(top:bot,left:right)+std(map(:)*randn(rr,cc));
        snapshot=0*map(top:bot,left:right);
    end
end

