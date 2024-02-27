function A010_get_spots(conditions,savename0)
exp_path='M:\tnw\bn\cd\Shared\Tisma\Bacillus subtilis chromosome imaging\Widefield_and_deconvolution\20230906_BSG5522_IPTG_timecourse\';
tif_path=[exp_path, '\tiffs\'];
jpeg_path=[exp_path, '\jpegs\'];
if ~isdir(jpeg_path), mkdir(jpeg_path); end

psf=3;
peeler_spotcount_wrapup = [];
peeler_spots_ori=[];
peeler_spots_ter=[];
oufti_spotcount_wrapup = [];
oufti_spots_ori=[];
oufti_spots_ter=[];
N_cond=length(conditions);
maxcells=5E10; %per frame
savename=[savename0,'_psf',num2str(psf)];

%1) work ROI

for ii=1:N_cond    
    mins=conditions(ii);
    minstr=num2str(mins);
    sourcename=['BSG5522_', minstr, 'min_c1_Stack'];
    spot_mat=[sourcename, '_SpotDetection_02.mat'];
    contour_mat=[sourcename, '_Corr_outlines.mat'];
    source=[exp_path,'\',spot_mat];
    load(source);
    cellL = cellList.meshData;
    cellcount = 0; % cell counter
    for frame = 1:length(cellL)
        %im_c1=(imread([tif_path, 'BSG5522_',minstr,'min_00',num2str(frame),'xy1c1.tif']));
        im_c2=(imread([tif_path, 'BSG5522_',minstr,'min_00',num2str(frame),'xy1c2.tif']));
        im_c3=(imread([tif_path, 'BSG5522_',minstr,'min_00',num2str(frame),'xy1c3.tif']));
        %im_c4=(imread([tif_path, 'BSG5522_',minstr,'min_00',num2str(frame),'xy1c4.tif']));
        Ncellstodo=min([ maxcells length(cellL{frame})]);
        for cell = 1:Ncellstodo
           if ~isempty (cellL{frame}{cell}) && isfield(cellL{frame}{cell},'spots') && isfield(cellL{frame}{cell},'spots2')...
                && ~isempty(cellL{frame}{cell}.spots.positions) && ~isempty(cellL{frame}{cell}.spots2.positions)... % usual conditions
                && length(cellL{frame}{cell}.spots.positions) >= 1 && length(cellL{frame}{cell}.spots2.positions) >= 1  % number of spots signal 1 and 2 must be more than 1  
                
                %% fetch oufti data
                N_ter=length(cellL{frame}{cell}.spots.x);  %ter number
                N_ori=length(cellL{frame}{cell}.spots2.x); %ori number                   
                box=cellL{frame}{cell}.box; %[bottom, left, width, height]
                roi_ter=im_c2(box(2):box(2)+box(4),box(1):box(1)+box(3));
                roi_ori=im_c3(box(2):box(2)+box(4),box(1):box(1)+box(3));
                ter_x_oufti=cellL{frame}{cell}.spots.x-box(1);
                ter_y_oufti=cellL{frame}{cell}.spots.y-box(2);
                ori_x_oufti=cellL{frame}{cell}.spots2.x-box(1);
                ori_y_oufti=cellL{frame}{cell}.spots2.y-box(2);
                cont_x_oufti=cellL{frame}{cell}.model(:,1)-box(1);
                cont_y_oufti=cellL{frame}{cell}.model(:,2)-box(2);
                
 
                if ~ isempty(cont_x_oufti)
                    cellcount = cellcount+1;  %workable cell, make mask
                    [mask,edge_mask,inner_mask]=make_mask_from_contour(roi_ter,cont_x_oufti,cont_y_oufti);
                    ter_med=median(roi_ter(edge_mask>0));
                    roi_ter=roi_ter-ter_med;
                    roi_ter(roi_ter<0)=0;
                    ori_med=median(roi_ori(edge_mask>0));
                    roi_ori=roi_ori-ori_med;
                    roi_ori(roi_ori<0)=0;
                    %some masking and smoothing work:
                    roi_ter=uint16(matrix_blur_jk(roi_ter.*mask,2));
                    roi_ori=uint16(matrix_blur_jk(roi_ori.*mask,2));
                    
                    %% oufti expand
                   
                    
                    spotprops_ori_oufti=expand_oufti_props(ori_x_oufti,ori_y_oufti,1*psf, roi_ori,mask,'ori'); 
                    [Nso_o,~]=size(spotprops_ori_oufti);
                    extrainfo_ori_oufti=repmat([mins frame cell],Nso_o,1);
                    oufti_spots_ori=[oufti_spots_ori ; [extrainfo_ori_oufti spotprops_ori_oufti]];
                    
                    spotprops_ter_oufti=expand_oufti_props(ter_x_oufti,ter_y_oufti,1*psf, roi_ter,mask,'ter');     
                    [Nso_t,~]=size(spotprops_ter_oufti);
                    extrainfo_ter_oufti=repmat([mins frame cell],Nso_t,1);
                    oufti_spots_ter=[oufti_spots_ter ; [extrainfo_ter_oufti spotprops_ter_oufti]];
                    oufti_spotcount_wrapup = [oufti_spotcount_wrapup ; [mins frame cell Nso_t Nso_o Nso_o/Nso_t]];
             
                    
                    %% peeler work 
                    [spotprops_ori_peeler,~]=matrix_peel_blobs_from_image(double(roi_ori),psf,1,2,'max_spots', 0.03,50,0);
                    [Nsp_o,~]=size(spotprops_ori_peeler);
                    extrainfo_ori=repmat([mins frame cell],Nsp_o,1);
                    peeler_spots_ori=[peeler_spots_ori ; [extrainfo_ori spotprops_ori_peeler]];
                    [spotprops_ter_peeler,~]=matrix_peel_blobs_from_image(double(roi_ter),psf,1,2,'max_spots', 0.03,50,0);
                    [Nsp_t,~]=size(spotprops_ter_peeler);
                    extrainfo_ter=repmat([mins frame cell],Nsp_t,1);
                    peeler_spots_ter=[peeler_spots_ter ; [extrainfo_ter spotprops_ter_peeler]];
                    %add summary of results to a single array:
                    peeler_spotcount_wrapup = [peeler_spotcount_wrapup ; [mins frame cell Nsp_t Nsp_o Nsp_o/Nsp_t]];
                    disp(cellcount);
                end
            end
        end
    end   
end

%wrapup 
%save spotdata:
    nme=[savename, '_spotdata.mat'];
    save(nme,'peeler_spots_ori', 'peeler_spots_ter', 'peeler_spotcount_wrapup', ...
             'oufti_spots_ori', 'oufti_spots_ter', 'oufti_spotcount_wrapup'); 

function spotprops=expand_oufti_props(xx,yy,psf, roi,mask,titl)
%aim for output: 
%spotprops=[idx pk x y Psf ThisSpotFraction CoveredFraction RelChange]];
roi=double(roi);
mask=double(mask);
Ls=length(xx);
%get local brightness
sum_I=sum(roi(:));
%[idx pk x y Psf ThisSpotFraction CoveredFraction RelChange]];
for ii=1:Ls
    count(ii)=ii;
    x=round(xx(ii));
    y=round(yy(ii));
    pk(ii)=roi(y,x);
    psfs(ii)=psf;
    bck_local(ii)=GetLocalBackgroundRing(roi,x,y,psf);
    temproi=roi-bck_local(ii);
    temproi(temproi<0)=0;
    gaussim_sum=sum(sum((pk(ii)*TwoDGaussNormPeak(temproi,x,y,psf))));
    localsum=GetLocalContentRing(temproi,x,y,psf);
    perc(ii)=gaussim_sum/sum_I;
end
%sorting
[perc,idx]=sort(perc, 'descend');
xx=xx(idx);
yy=yy(idx);
pk=pk(idx);
psfs=psfs(idx);
cumperc=cumsum(perc);
%spotprops=[idx pk x y Psf ThisSpotFraction CoveredFraction RelChange]];
spotprops=[count' pk' xx' yy' psfs' perc' cumperc' NaN*count' bck_local'];
if 0
    figure(93);
    spotmaxscale=120;
    pcolor(roi.*mask); shading flat; axis equal; axis tight; colormap bone;  hold on;
    plot(xx+0.5, yy+0.5,'bo', 'MarkerSize', 8, 'LineWidth',1);
    title(titl);
    for ii=1:Ls
        x=xx(ii);
        y=yy(ii);
        p=perc(ii)
        mrksize=p*spotmaxscale/100;
        plot(x,y,'yo','MarkerSize',max([mrksize 1]),'LineWidth',1);
        text(x+10,y,num2str(round(p)), 'Color', 'y');
    end
    hold off;
    pause(0.1);
    [~]=ginput(1);
end





