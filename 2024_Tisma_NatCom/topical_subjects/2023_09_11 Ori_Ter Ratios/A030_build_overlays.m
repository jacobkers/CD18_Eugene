function A030_build_overlays(conditions,savename0)
exp_path='M:\tnw\bn\cd\Shared\Tisma\Bacillus subtilis chromosome imaging\Widefield_and_deconvolution\20230906_BSG5522_IPTG_timecourse\';
tif_path=[exp_path, '\tiffs\'];
jpeg_path=[exp_path, '\jpegs\'];
if ~isdir(jpeg_path), mkdir(jpeg_path); end

psf=3;
treshold_amplifier=1.0;
treshold_wrapup_ori=[];
treshold_wrapup_ter=[];
N_cond=length(conditions);
skips_jpeg=1;  %for jpgs
maxcells=3E10;     %per frame
savename=[savename0,'_psf',num2str(psf)];

nme=[savename, '_spotdata.mat'];
load(nme,'all_spots_ori', 'all_spots_ter', 'spotcount_wrapup', 'treshold_wrapup_ori', 'treshold_wrapup_ter'); 


%1) load per cell
for ii=1:N_cond
    mins=conditions(ii);
    minstr=num2str(mins);
    itx=find(treshold_wrapup_ori(:,1)==mins);
    treshold_ori=treshold_amplifier*treshold_wrapup_ori(itx,5);
    treshold_ter=treshold_amplifier*treshold_wrapup_ter(itx,5);
    
    spots_ori_thiscondition=all_spots_ori(all_spots_ori(:,1)==mins,:);
    spots_ter_thiscondition=all_spots_ter(all_spots_ter(:,1)==mins,:);   
    spot_mat=['BSG5522_', minstr, 'min_c1_Stack_SpotDetection_02.mat'];
    contour_mat=['BSG5522_', minstr, 'min_c1_Stack_Corr_outlines.mat'];
    source=[exp_path,'\',spot_mat];
    load(source);
    cellL = cellList.meshData;
    cellcount = 0; % cell counter
    for frame = 1:length(cellL)
        %BSG5522_0min_002xy1c2
        im_c1=(imread([tif_path, 'BSG5522_',minstr,'min_00',num2str(frame),'xy1c1.tif']));
        im_c2=(imread([tif_path, 'BSG5522_',minstr,'min_00',num2str(frame),'xy1c2.tif']));
        im_c3=(imread([tif_path, 'BSG5522_',minstr,'min_00',num2str(frame),'xy1c3.tif']));
        im_c4=(imread([tif_path, 'BSG5522_',minstr,'min_00',num2str(frame),'xy1c4.tif']));
        Ncellstodo=min([ maxcells length(cellL{frame})]);
        for cell = 1:Ncellstodo
               if   ~isempty (cellL{frame}{cell}) && isfield(cellL{frame}{cell},'spots') && isfield(cellL{frame}{cell},'spots2')...
                    && ~isempty(cellL{frame}{cell}.spots.positions) && ~isempty(cellL{frame}{cell}.spots2.positions)... % usual conditions
                    && length(cellL{frame}{cell}.spots.positions) >= 1 && length(cellL{frame}{cell}.spots2.positions) >= 1  % number of spots signal 1 and 2 must be more than 1
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
                            cellcount = cellcount+1;  %workable cell
                            cell_idx=find(spots_ori_thiscondition(:,1)==mins&...
                                          spots_ori_thiscondition(:,2)==frame&...
                                          spots_ori_thiscondition(:,3)==cell);
                            spotprops_thiscell_ori=spots_ori_thiscondition(cell_idx,:);
                            spotpropsthiscell_ter=spots_ter_thiscondition(cell_idx,:);                      
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
                            
                             %show:
                            if mod(cellcount, skips_jpeg) == 0
                            figure(79);
                            subplot(1,2,1);
                                pcolor(roi_ter.*mask); shading flat; axis equal; axis tight; colormap bone;  hold on;
                                plot(ter_x_oufti+0.5, ter_y_oufti+0.5,'bo', 'MarkerSize', 8, 'LineWidth',1);
                                plot(cont_x_oufti+0.5,cont_y_oufti+0.5,'m-');
                                plot_spots(spotpropsthiscell_ter,treshold_ter);
                                %plot(ori_x-box(1)+0.5, ori_y-box(2)+0.5,'ro','MarkerSize', 6);
                                title('ter channel');
                                hold off;
                            subplot(1,2,2);
                                pcolor(roi_ori.*mask); shading flat; axis equal; axis tight; colormap bone;  hold on;
                                %plot(ter_x-box(1)+0.5, ter_y-box(2)+0.5,'bo', 'MarkerSize', 6);
                                plot(ori_x_oufti+0.5, ori_y_oufti+0.5,'ro','MarkerSize', 8, 'LineWidth',1);
                                plot_spots(spotprops_thiscell_ori,treshold_ori);
                                plot(cont_x_oufti+0.5,cont_y_oufti+0.5,'w-');
                                title('ori channel');
                                hold off;;
                                pause(0.1);
                            nme=['time',minstr,'mins_frame', num2str(frame,'%03.0f'), 'cell',num2str(cell,'%04.0f'),'.jpg'];
                            saveas(gcf,[jpeg_path nme]);
                            pause(0.1);
                            %[~]=ginput(1);
                            end
                  end
              end
        end
    end    
end

    function plot_spots(spotprops, treshold);
        spotmaxscale=50;  %largest spot
        sel=find(spotprops(:,9)>treshold);
        if ~isempty(sel)
            bright_spotprops=spotprops(sel,:);
            Nmax=15;
            [LS,~]=size(bright_spotprops);
            Nspots=min([Nmax LS]);
            for sp=1:Nspots
                mrksize=ceil(bright_spotprops(sp,9)*spotmaxscale);
                x=bright_spotprops(sp,6);
                y=bright_spotprops(sp,7);
                plot(x,y,'yo','MarkerSize',max([mrksize 1]),'LineWidth',1);
                hold on; 
            end
        end

