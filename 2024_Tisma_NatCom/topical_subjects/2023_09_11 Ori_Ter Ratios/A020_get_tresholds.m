function A020_get_tresholds(conditions,savename0)
exp_path='M:\tnw\bn\cd\Shared\Tisma\Bacillus subtilis chromosome imaging\Widefield_and_deconvolution\20230906_BSG5522_IPTG_timecourse\';
tif_path=[exp_path, '\tiffs\'];
jpeg_path=[exp_path, '\jpegs\'];
if ~isdir(jpeg_path), mkdir(jpeg_path); end

psf=3;

treshold_wrapup_ori_peel=[];
treshold_wrapup_ter_peel=[];
treshold_wrapup_ori_oufti=[];
treshold_wrapup_ter_oufti=[];
N_cond=length(conditions);

savename=[savename0,'_psf',num2str(psf)];

nme=[savename, '_spotdata.mat'];
load(nme,'peeler_spots_ori', 'peeler_spots_ter', ...
             'oufti_spots_ori', 'oufti_spots_ter'); 


%% 1) work conditions to find typical tresholds per conditions
for ii=1:N_cond
    mins=conditions(ii);
    %peeler
    spots_ori_peeler_thiscondition=peeler_spots_ori(peeler_spots_ori(:,1)==mins,:);
    spots_ter_peeler_thiscondition=peeler_spots_ter(peeler_spots_ter(:,1)==mins,:);
    %find treshold ori using all detected peeler spots:
    figure(168);
    [N_typical_ori_peel, Perc_limit_ori_peel,cellcount_peel]= get_treshold_info(spots_ori_peeler_thiscondition,1);
    title('ori_peeler');
    figure(169);
    [N_typical_ter_peel, Perc_limit_ter_peel,cellcount_peel]= get_treshold_info(spots_ter_peeler_thiscondition,1);
    title('ter_peeler');
    treshold_wrapup_ori_peel=[treshold_wrapup_ori_peel; ...
        [mins psf cellcount_peel N_typical_ori_peel Perc_limit_ori_peel]];
    treshold_wrapup_ter_peel=[treshold_wrapup_ter_peel; ...
        [mins psf cellcount_peel N_typical_ter_peel Perc_limit_ter_peel]];  
    
    %oufti
    spots_ori_oufti_thiscondition=oufti_spots_ori(oufti_spots_ori(:,1)==mins,:);
    spots_ter_oufti_thiscondition=oufti_spots_ter(oufti_spots_ter(:,1)==mins,:);
    %find treshold ori using all detected peeler spots:
    figure(180);
    [N_typical_ori_oufti, Perc_limit_ori_oufti,cellcount_oufti]= get_treshold_info(spots_ori_oufti_thiscondition,1);
    title('ori_oufti');
    figure(181);
    [N_typical_ter_oufti, Perc_limit_ter_oufti,cellcount_oufti]= get_treshold_info(spots_ter_oufti_thiscondition,1);
    title('ter_oufti');
    
     treshold_wrapup_ori_oufti=[treshold_wrapup_ori_oufti; ...
        [mins psf cellcount_oufti N_typical_ori_oufti Perc_limit_ori_oufti]];
     treshold_wrapup_ter_oufti=[treshold_wrapup_ter_oufti; ...
        [mins psf cellcount_oufti N_typical_ter_oufti Perc_limit_ter_oufti]]; 

end


%% wrapup : save new data
%save spotdata:
    nme=[savename, '_spotdata.mat'];
    save(nme, 'treshold_wrapup_ori_peel', 'treshold_wrapup_ter_peel',...
              'treshold_wrapup_ori_oufti', 'treshold_wrapup_ter_oufti', '-append'); 

%saving processded data: excel 
excelname=[savename, '.xlsx'];
ColNames=[{'minutes'} , {'psf'}, {'cellcount'} ,{'N_percell_limit'}, {'Perc_limit'}];
xlswrite(excelname,ColNames,'ori','A1');
xlswrite(excelname,treshold_wrapup_ori_peel,'ori','A2');
xlswrite(excelname,ColNames,'ter','A1');
xlswrite(excelname,treshold_wrapup_ter_peel,'ter','A2');

function [N_typical, Perc_limit,cellcount]= get_treshold_info(all_spots,show)
    %mins frame cell
    cellcount=length(find(diff(all_spots(:,3))>0));
    peakz_asc=sort(all_spots(:,9),'ascend');
    peakz_des=sort(all_spots(:,9),'descend');
    [tres_perc,xp]=JKD1_PRF_Find_treshold_MD_1D(peakz_asc,0);
    
    axz=(1:length(peakz_asc))'/cellcount;
    idx_des=(length(peakz_asc)-xp+1);
    N_typical=idx_des/cellcount;
    Perc_limit=tres_perc;
    if show
        %figure;
        plot(axz,100*peakz_des, '-'); hold on;
        plot(axz(idx_des),100*peakz_des(idx_des), 'o');
        ylabel('spot %')
        xlabel('spotindex/Ncells')
        pause(0.5);
        [~]=ginput(1);
    end
    
    function plot_spots(spotprops_ori);
        Nmax=7;
        spotmaxscale=24;
        scale_ori=24/max(spotprops_ori(:,8));
        [LS,~]=size(spotprops_ori);
        Nspots=min([Nmax LS])
        for sp=1:Nspots
            mrksize=ceil(spotprops_ori(sp,8)*spotmaxscale);
            x=spotprops_ori(sp,5);
            y=spotprops_ori(sp,6);
            plot(x,y,'yo','MarkerSize',max([mrksize 1]),'LineWidth',1);
            hold on; 
        end

