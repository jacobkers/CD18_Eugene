function A026_get_histograms_oufti(conditions,savename0)
close all;
exp_path='M:\tnw\bn\cd\Shared\Tisma\Bacillus subtilis chromosome imaging\Widefield_and_deconvolution\20230906_BSG5522_IPTG_timecourse\';
tif_path=[exp_path, '\tiffs\'];
jpeg_path=[exp_path, '\jpegs\'];
if ~isdir(jpeg_path), mkdir(jpeg_path); end


all_spots_ori=[];
all_spots_ter=[];
treshold_wrapup_ori=[];
treshold_wrapup_ter=[];
N_cond=length(conditions);
treshold_amplifier=1.0;
maxcells=5E6; %per frame
savename=['A026_', savename0];

nme=[savename, '_spotdata.mat'];
load(['oufti_data\',nme], 'oufti_spots_ori',  'oufti_spots_ter',  'treshold_wrapup_ori_oufti'); 

goodspots_ori=[];
goodspots_ter=[];
good_spotcount_wrapup=[];

%set up output data:
binz_oriter=0:1:13;
binz_ratio=0:1:10;
barplotdata_ori=zeros(length(binz_oriter),N_cond+1);
barplotdata_ter=zeros(length(binz_oriter),N_cond+1);
barplotdata_ratio=zeros(length(binz_ratio),N_cond+1);
headers_oriter=[{'number of foci'}];
headers_ratio=[{'ori_ter ratio'}];
for ci=1:N_cond
    headers_oriter{ci+1}=num2str(conditions(ci));
    headers_ratio{ci+1}=num2str(conditions(ci));
end
barplotdata_ori(:,1)=binz_oriter;
barplotdata_ter(:,1)=binz_oriter;
barplotdata_ratio(:,1)=binz_ratio;

%% evaluate per workable cell:
for ci=1:N_cond
    mins=conditions(ci);
    minstr=num2str(mins);
    itx=find(treshold_wrapup_ori_oufti(:,1)==mins);
    %tresholds:
    treshold_ori=0.15;  
    treshold_ter=0.15;
    %all spots
    spots_ori_thiscondition=oufti_spots_ori(oufti_spots_ori(:,1)==mins,:);
    spots_ter_thiscondition=oufti_spots_ter(oufti_spots_ter(:,1)==mins,:);   
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

    %build panels for this condition
    spotcount_wrapup_thiscondition=good_spotcount_wrapup(good_spotcount_wrapup(:,1)==mins,:);
    N_spot_ter(ci)=mean(spotcount_wrapup_thiscondition(:,4));
    N_spot_ori(ci)=mean(spotcount_wrapup_thiscondition(:,5));
    N_spot_ratio(ci)=mean(spotcount_wrapup_thiscondition(:,6));
    
       
    figure(112);
    Hist_N_ter = hist(spotcount_wrapup_thiscondition(:,4),binz_oriter);
    Hist_N_ori = hist(spotcount_wrapup_thiscondition(:,5),binz_oriter);
    Hist_N_ori_ter_ratio = hist(good_spotcount_wrapup(:,6),binz_ratio);
    
    barplotdata_ori(:,ci+1)=Hist_N_ori';
    barplotdata_ter(:,ci+1)=Hist_N_ter';
    barplotdata_ratio(:,ci+1)=Hist_N_ori_ter_ratio';
    
    %figure;
    subplot(N_cond,3,1+3*(ci-1)); bar(binz_oriter,Hist_N_ori, 'r');
    legend('ori')
    title(minstr)
    subplot(N_cond,3,2+3*(ci-1)); bar(binz_oriter,Hist_N_ter, 'b');
    legend('ter')
    subplot(N_cond,3,3+3*(ci-1)); bar(binz_ratio,Hist_N_ori_ter_ratio, 'w');
    legend('ratio')
    dum=1;
end

%saving processded data: jpg
figure(112);
jpgname1=[savename, '_histograms.jpeg'];
saveas(gcf,jpgname1);


figure(105);
subplot(1,2,1);
plot(conditions,N_spot_ori,'ro-'); hold on;
plot(conditions,N_spot_ter, 'bo-'); hold on;
legend('ori','ter');
subplot(1,2,2);
plot(conditions,N_spot_ratio,'o-'); hold on;
legend('ratio');


jpgname2=[savename, '_counts_and_ratios.jpeg'];
saveas(gcf,jpgname2);


%% build exel output
%build excel table of barplots
xlswrite([savename, '_histograms.xls'], headers_oriter, 'ori', 'A1');
xlswrite([savename, '_histograms.xls'], barplotdata_ori, 'ori', 'A2');
xlswrite([savename, '_histograms.xls'], headers_oriter, 'ter', 'A1');
xlswrite([savename, '_histograms.xls'], barplotdata_ter, 'ter', 'A2');
xlswrite([savename, '_histograms.xls'], headers_oriter, 'ratio', 'A1');
xlswrite([savename, '_histograms.xls'], barplotdata_ratio, 'ratio', 'A2');

%build excel table of barplots
data_wrapup1=[conditions',N_spot_ori',N_spot_ter',N_spot_ratio']; 
headers_wrapup1=[{'time, mins'},{'#ori'},{'#ter'}, {'ratio'}];
xlswrite([savename, '_counts_and_ratios.xls'], headers_wrapup1, 'counts_and_ratios', 'A1');
xlswrite([savename, '_counts_and_ratios.xls'], data_wrapup1, 'counts_and_ratios', 'A2');

dum=1;







