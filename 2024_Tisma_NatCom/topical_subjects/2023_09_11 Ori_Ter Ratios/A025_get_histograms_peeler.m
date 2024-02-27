function A025_get_histograms(conditions,savename0)
close all;
exp_path='M:\tnw\bn\cd\Shared\Tisma\Bacillus subtilis chromosome imaging\Widefield_and_deconvolution\20230906_BSG5522_IPTG_timecourse\';
tif_path=[exp_path, '\tiffs\'];
jpeg_path=[exp_path, '\jpegs\'];
if ~isdir(jpeg_path), mkdir(jpeg_path); end

psf=3;
all_spots_ori=[];
all_spots_ter=[];
treshold_wrapup_ori=[];
treshold_wrapup_ter=[];
N_cond=length(conditions);
treshold_amplifier=1.0;
maxcells=5E6; %per frame
savename=[savename0,'_psf',num2str(psf)];

nme=[savename, '_spotdata.mat'];
load(nme, 'peeler_spots_ori', 'peeler_spots_ter', 'treshold_wrapup_ori_peel', 'treshold_wrapup_ter_peel',...
           'oufti_spots_ori', 'oufti_spots_ter','treshold_wrapup_ori_oufti', 'treshold_wrapup_ter_oufti'); 



goodspots_ori=[];
goodspots_ter=[];
good_spotcount_wrapup=[];

%% evaluate per workable cell:
for ii=1:N_cond
    mins=conditions(ii);
    minstr=num2str(mins);
    itx=find(treshold_wrapup_ori_peel(:,1)==mins);
    %tresholds:
    treshold_ori=treshold_amplifier*treshold_wrapup_ori_peel(itx,5);
    treshold_ter=treshold_amplifier*treshold_wrapup_ter_peel(itx,5);
    %all spots
    spots_ori_thiscondition=peeler_spots_ori(peeler_spots_ori(:,1)==mins,:);
    spots_ter_thiscondition=peeler_spots_ter(peeler_spots_ter(:,1)==mins,:);   
    spot_mat=['BSG5522_', minstr, 'min_c1_Stack_SpotDetection_02.mat'];

    source=[exp_path,'\',spot_mat];
    load(source);
    cellL = cellList.meshData;
    cellcount = 0; % cell counter
    for frame = 1:length(cellL)      
        Ncellstodo=min([ maxcells length(cellL{frame})]);
        for cell = 1:Ncellstodo
               if   ~isempty (cellL{frame}{cell}) && isfield(cellL{frame}{cell},'spots') && isfield(cellL{frame}{cell},'spots2')...
                    && ~isempty(cellL{frame}{cell}.spots.positions) && ~isempty(cellL{frame}{cell}.spots2.positions)... % usual conditions
                    && length(cellL{frame}{cell}.spots.positions) >= 1 && length(cellL{frame}{cell}.spots2.positions) >= 1  % number of spots signal 1 and 2 must be more than 1
                    box=cellL{frame}{cell}.box; %[bottom, left, width, height]
                    cont_x_oufti=cellL{frame}{cell}.model(:,1)-box(1);
                    cont_y_oufti=cellL{frame}{cell}.model(:,2)-box(2);    
                    if ~ isempty(cont_x_oufti)
                            cellcount = cellcount+1;  %workable cell
                            cell_idx_ori=find(spots_ori_thiscondition(:,1)==mins&...
                                          spots_ori_thiscondition(:,2)==frame&...
                                          spots_ori_thiscondition(:,3)==cell);
                            cell_idx_ter=find(spots_ter_thiscondition(:,1)==mins&...
                                          spots_ter_thiscondition(:,2)==frame&...
                                          spots_ter_thiscondition(:,3)==cell);
                            spotprops_thiscell_ori=spots_ori_thiscondition(cell_idx_ori,:);
                            spotprops_thiscell_ter=spots_ter_thiscondition(cell_idx_ter,:); 
                            %[cond frame cell spotcount Peak Xpos Ypos Psf Fraction CoveredFraction RelChange]];
                            good_ter_idx= find(spotprops_thiscell_ter(:,9)>treshold_ter);
                            good_ori_idx= find(spotprops_thiscell_ori(:,9)>treshold_ori);
                            if ~isempty(good_ter_idx)&~isempty(good_ori_idx)
                                goodspots_ori_thiscell=spotprops_thiscell_ori(good_ori_idx);
                                goodspots_ter_thiscell=spotprops_thiscell_ter(good_ter_idx);
                                goodspots_ori=[goodspots_ori ; goodspots_ori_thiscell ];
                                goodspots_ter=[goodspots_ter ; goodspots_ter_thiscell ];
                                [Nsp_o,~]=size(goodspots_ori_thiscell);
                                [Nsp_t,~]=size(goodspots_ter_thiscell);
                                good_spotcount_wrapup = [good_spotcount_wrapup ; [mins frame cell Nsp_t Nsp_o Nsp_o/Nsp_t]];
                            end
                            
                  end
              end
        end
    end 
    
    dum=1
    
    %build panels for this condition
    spotcount_wrapup_thiscondition=good_spotcount_wrapup(good_spotcount_wrapup(:,1)==mins,:);
    N_spot_ter(ii)=mean(spotcount_wrapup_thiscondition(:,4));
    N_spot_ori(ii)=mean(spotcount_wrapup_thiscondition(:,5));
    N_spot_ratio(ii)=mean(spotcount_wrapup_thiscondition(:,6));
    figure(112);
    binz1=0:1:13;
    binz2=0:1:10;
    Hist_N_ter = hist(spotcount_wrapup_thiscondition(:,4),binz1);
    Hist_N_ori = hist(spotcount_wrapup_thiscondition(:,5),binz1);
    Hist_N_ori_ter_ratio = hist(good_spotcount_wrapup(:,6),binz2);
    %figure;
    subplot(N_cond,3,1+3*(ii-1)); bar(binz1,Hist_N_ori, 'r');
    legend('ori')
    title(minstr)
    subplot(N_cond,3,2+3*(ii-1)); bar(binz1,Hist_N_ter, 'b');
    legend('ter')
    subplot(N_cond,3,3+3*(ii-1)); bar(binz2,Hist_N_ori_ter_ratio, 'w');
    legend('ratio')
    dum=1;
end

%saving processded data: jpg
figure(112);
jpgname1=[savename, '_histograms.jpeg'];
%saveas(gcf,jpgname1);

figure(105);
subplot(1,2,1);
plot(conditions,N_spot_ori,'ro-'); hold on;
plot(conditions,N_spot_ter, 'bo-'); hold on;
legend('ori','ter');
subplot(1,2,2);
plot(conditions,N_spot_ratio,'o-'); hold on;
legend('ratio');


jpgname2=[savename, '_ounts_and_ratios.jpeg'];
saveas(gcf,jpgname2);






