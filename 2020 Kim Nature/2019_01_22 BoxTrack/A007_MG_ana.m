%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script analyzes kymograph generated from MG_gen
% Originally written by Marijn van Loenhout;
% edited/adapted by Jacob Kers 2013
% edited/adapted by SHKIM DEC2013
% edited/adapted by SHKIM OCT2015
% Version 8 by SHKIM on JUL2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------

%% Init Parameters
%     clear all;
warning off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analysis of Two KymoGraphs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% plot the kymo of coiled DNA
[Npixels,tmax]=size(KymoBgBlCoil);
figure(31);clf;
subplot(3,4,[5 6 9 10]);
imagesc(KymoBgBlCoil');
caxis([0 5]);   % 5 is arbitrary number. This gives the same color scale in intensity kymograph
title('Intensity kymo');
xlabel('Position (pixel)');
ylabel('Time (frame)');
colorbar;



%% calculate mean intensity profile of supercoild and relaxed molecule-------------
plecI=mean(KymoBgBlCoil,2);  % normalzed during bleach correction
plecI_sigma=std(KymoBgBlCoil,0,2);  % normalzed during bleach correction
nickedI=mean(KymoBgBlNick,2);
nickedI_sigma=std(KymoBgBlNick,0,2);
x_px=0:length(nickedI)-1;                %axis nicked in original pixels

%% make interpolated DNA intensity profile
x_ip=(0:1/initval.interpolation_ratio:length(nickedI)-1);    %same axis interpolated in 100x as much equal parts
nickedI_ip=interp1(x_px,nickedI,x_ip,'spline');     %nicked profile re-mapped on new axisv
plecI_ip=interp1(x_px,plecI,x_ip,'spline');     %nicked profile re-mapped on new axisv


%% plot average intensity profile
figure(31); 
subplot(3,4,[1 2 3 4]);
errorbar(x_px,nickedI,nickedI_sigma,'b.');hold on
fhdle1=plot(x_ip,nickedI_ip,'b-');
title('Averaged intensity profile');
ylabel('Intensity');
errorbar(x_px,plecI,plecI_sigma,'r.');
fhdle2=plot(x_ip,plecI_ip,'r-');
axis tight;
xlabel('Position (pixel)');
legend([fhdle1, fhdle2],'relaxed','coiled','Location','best'); legend('boxoff');

hold off;





%% determine DNA attachment; get length of tether---------------------------

% if find_end_manual==2
%     thresh_DNAend=0.5;
% 
%     nickedI_tmp=median(KymoBgBlNick,2);
%     plecI_tmp=median(KymoBgBlCoil,2);
%     DNAinfo.tresh0=thresh_DNAend*mean(nickedI_tmp(nickedI_tmp>median(nickedI_tmp)/2));
%     DNAinfo.tresh=thresh_DNAend*mean(plecI_tmp(plecI_tmp>median(plecI_tmp)/2));
% 
%     DNA_start_Nick_px_ip=find(nickedI_ip>DNAinfo.tresh0,1);
%     DNA_start_Coil_avr_ip=find(plecI_ip>DNAinfo.tresh,1);
%     DNA_end_Nick_px_ip=find(nickedI_ip>DNAinfo.tresh0,1,'last');
%     DNA_end_Coil_avr_ip=find(plecI_ip>DNAinfo.tresh,1,'last');
% 
%     DNAinfo.start_Nick_px=DNA_start_Nick_px_ip/initval.interpolation_ratio;
%     DNAinfo.start_Coil_avr=DNA_start_Coil_avr_ip/initval.interpolation_ratio;
%     DNAinfo.end_Nick_px=DNA_end_Nick_px_ip/initval.interpolation_ratio;
%     DNAinfo.end_Coil_avr=DNA_end_Coil_avr_ip/initval.interpolation_ratio;
% else
% 	DNA_start_Nick_px_ip=DNAinfo.start_Nick_px*initval.interpolation_ratio;
%     DNA_start_Coil_avr_ip=DNAinfo.start_Coil_avr*initval.interpolation_ratio;
%     DNA_end_Nick_px_ip=DNAinfo.end_Nick_px*initval.interpolation_ratio;
%     DNA_end_Coil_avr_ip=DNAinfo.end_Coil_avr*initval.interpolation_ratio;
% end
% 
% subplot(3,4,[1 2 3 4]);
% [bar_size,~]=max(nickedI);
% bar_size=bar_size*2;
% line([DNAinfo.end_Nick_px,DNAinfo.end_Nick_px],[0,bar_size],'color','b');
% line([DNAinfo.start_Nick_px,DNAinfo.start_Nick_px],[0,bar_size],'color','b');
% line([0,Npixels],[DNAinfo.tresh0,DNAinfo.tresh0],'color','b');
% line([DNAinfo.end_Coil_avr,DNAinfo.end_Coil_avr],[0,bar_size],'color','r');
% line([DNAinfo.start_Coil_avr,DNAinfo.start_Coil_avr],[0,bar_size],'color','r');
% line([0,Npixels],[DNAinfo.tresh,DNAinfo.tresh],'color','r');
% hold off



%%


% find_end_manual=menu('End determination','manual','automatic');

% if find_end_manual==1
    figure(31)
    subplot(3,4,[1 2 3 4]);
    errorbar(x_px,nickedI,nickedI_sigma,'b.-');
    title('Averaged intensity profile');
    ylabel('Intensity');
    ylim([0 3])
    xlim([5 max(x_px)-5])
    hold on;
    errorbar(x_px,plecI,plecI_sigma,'r.-');
%     axis tight;
    xlabel('Position (pixel)');
    legend('relaxed','colied','Location','best'); legend('boxoff');
    
   %% 
    % nick DNA
    title('Click for relaxed DNA');
    [ginpx,~]=ginput(2);
    
    DNAinfo.start_Nick_px_ip=ginpx(1)*initval.interpolation_ratio;
    DNAinfo.end_Nick_px_ip=ginpx(2)*initval.interpolation_ratio;
    DNAinfo.start_Nick_px=round(ginpx(1))+1;
    DNAinfo.end_Nick_px=round(ginpx(2))+1;
    
    % coil DNA
    title('Click for coiled DNA');
    [ginpx,~]=ginput(2);
    DNAinfo.start_Coil_avr_ip=ginpx(1)*initval.interpolation_ratio;
    DNAinfo.end_Coil_avr_ip=ginpx(2)*initval.interpolation_ratio;
    DNAinfo.start_Coil_avr=round(ginpx(1))+1;
    DNAinfo.end_Coil_avr=round(ginpx(2))+1;
    
    if DNAinfo.start_Nick_px<1, DNAinfo.start_Nick_px=1;    end
    if DNAinfo.start_Coil_avr<1, DNAinfo.start_Coil_avr=1;    end
    if DNAinfo.end_Nick_px>Npixels, DNAinfo.end_Nick_px=Npixels;    end
    if DNAinfo.end_Coil_avr>Npixels, DNAinfo.end_Coil_avr=Npixels;    end
    
%     DNAinfo.start_Nick_px=DNAinfo.start_Nick_px_ip/initval.interpolation_ratio;
%     DNAinfo.start_Coil_avr=DNAinfo.start_Coil_avr_ip/initval.interpolation_ratio;
%     DNAinfo.end_Nick_px=DNAinfo.end_Nick_px_ip/initval.interpolation_ratio;
%     DNAinfo.end_Coil_avr=DNAinfo.end_Coil_avr_ip/initval.interpolation_ratio;
% end



DNAinfo.len_Nick_px=DNAinfo.end_Nick_px-DNAinfo.start_Nick_px;
DNAinfo.len_Nick_um=DNAinfo.len_Nick_px*(initval.Px2um);
DNAinfo.len_Coil_px=DNAinfo.end_Coil_avr-DNAinfo.start_Coil_avr;
DNAinfo.len_Coil_um=DNAinfo.len_Coil_px*(initval.Px2um);
DNAinfo.superFraction=(DNAinfo.len_Nick_px-DNAinfo.len_Coil_px)/DNAinfo.len_Nick_px;
DNAinfo.N_base_per_px_Nick=DNAinfo.DNAlen_bp/DNAinfo.len_Nick_px/1000;   % Unit:kb



disp(['DNA length Nicked: ' num2str(DNAinfo.len_Nick_px) ' px (' num2str(DNAinfo.len_Nick_um) ' um)']);
disp(['DNA length Coiled: ' num2str(DNAinfo.len_Coil_px) ' px (' num2str(DNAinfo.len_Coil_um) ' um)']);




%--------------------------------------------------------------------------
%% calculate cummulative intensity profile of relaxed DNA
% nickedCum_ip=cumsum(nickedI_ip(DNA_start_Nick_px_ip:DNA_end_Nick_px_ip));
nickedCum_ip=cumsum(nickedI_ip(DNAinfo.start_Nick_px_ip:DNAinfo.end_Nick_px_ip));

% nickedCum_ip=cumsum(nickedI_ip);
nickedCum_ip=nickedCum_ip-min(nickedCum_ip);
nickedCum_ip=(nickedCum_ip)/max(nickedCum_ip);

x_px_cut=1:(DNAinfo.end_Nick_px-DNAinfo.start_Nick_px+1);
x_ip_cut=1:1/initval.interpolation_ratio:length(nickedCum_ip);    %same axis interpolated in 100x as much equal parts

%% ------------------------------------------------------------------------
% map the cummulative intensity profile of coiled molecule onto relaxed ---
%     disp('mapping the relaxed to the coiled DNA');
Kymo_DNAdensity=[];
Kymo_pxdensity=[];
Kymo_angle_ratio=[];
% DNAinfo.start_coil_px=(1:tmax)*0;
% DNAinfo.end_coil_px=(1:tmax)*0;
% DNAinfo.len_coil_px=(1:tmax)*0;

% tmax=500;
for fr_id=1:tmax
    %% get intensity profile of the current frame
    int_coil = KymoBgBlCoil(:,fr_id)';
    int_coil = smooth(int_coil,3,'sgolay'); % for DNA density analysis    
%     int_coil_ip=interp1(x_px,int_coil,x_ip,'spline');     %nicked profile re-mapped on new axisv
       
    %% build normalized cum int profile
    CumPlecI=cumsum(int_coil(ceil(DNAinfo.start_Coil_avr):floor(DNAinfo.end_Coil_avr)));
    CumPlecI=CumPlecI-min(CumPlecI);
    CumPlecI=(CumPlecI)/max(CumPlecI);
% 
%     figure(3255);
%     subplot(1,2,1);
%     plot(x_px,nickedI,'b.');hold on
%     fhdle1=plot(x_ip,nickedI_ip,'b-');
%     title('Averaged intensity profile');
%     ylabel('Intensity');
%     plot(x_px,plecI,'r.');
%     fhdle2=plot(x_ip,plecI_ip,'r-');
%     axis tight;
%     xlabel('Position (pixel)');
%     legend([fhdle1, fhdle2],'relaxed','coiled','Location','best'); legend('boxoff');
%     hold off;
%     
%     subplot(1,2,2);
%     plot((1:length(nickedCum_ip))/initval.interpolation_ratio,nickedCum_ip,'b');    hold on;
%     plot((1:length(CumPlecI))-1,CumPlecI,'r');  hold off;
%     axis tight;
%     title(['frame #' num2str(fr_id) ' / ' num2str(tmax)]);
%     drawnow
    
%         %% Map onto the relaxed DNA
        neo_px_id=1;
        id_map=zeros(1,length(int_coil));
%         id_map=KymoBgBlCoil(:,fr_id)';
%         if find_end_manual==1
%             for px_id=DNAinfo.start_Nick_px+1:DNAinfo.end_Nick_px-1
            for px_id=DNAinfo.start_Coil_avr+1:DNAinfo.end_Coil_avr-1
                neo_px_id=neo_px_id+1;
                [~,nick10ki]=min((nickedCum_ip-CumPlecI(neo_px_id)).^2); %look for best-matching summed value on nicked curve
                id_map(px_id)=nick10ki/initval.interpolation_ratio;
            end
%         else
%             for px_id=ceil(DNAinfo.start_Coil_avr)+1:floor(DNAinfo.end_Coil_avr)-1
%                 neo_px_id=neo_px_id+1;
%                 [~,nick10ki]=min((nickedCum_ip-CumPlecI(neo_px_id)).^2); %look for best-matching summed value on nicked curve
%                 id_map(px_id)=nick10ki/initval.interpolation_ratio;  %note: 100x precision
%             end
%         end
        index_jump=id_map-[0 id_map(1:end-1)];  %derivative curve: gives 'index jumps'
        index_jump=index_jump*DNAinfo.N_base_per_px_Nick;   % convert index jump (in px space) to # base (basepair space)
        index_jump(index_jump<0)=0; % to remove the final big jump due to the length difference of the two DNA
        
        
        Kymo_DNAdensity=[Kymo_DNAdensity;index_jump];
%         Kymo_DNAdensity=[Kymo_DNAdensity;id_map];
    
end
% Kymo_DNAdensity=KymoBgBlCoil';

% plot the DNA density kymograph
figure(31);
subplot(3,4,[7 8 11 12]);
imagesc(Kymo_DNAdensity);
caxis([0 initval.pt_density_imscale]);
xlabel('Position (pixel)');
title('DNA density map (unit: kb)');
colorbar;



%% save DNA density kymograph
cur_time=clock;
cur_time(6)=floor(cur_time(6));
kymo_output_dir=['kymo_' num2str(cur_time(1)) num2str(cur_time(2)) num2str(cur_time(3)) ...
    '_' num2str(cur_time(4)) 'h' num2str(cur_time(5)) 'm' num2str(cur_time(6),2) 's'];
if MG_auto_run
    kymo_output_dir=['re3_' kymo_output_dir];
end
    
mkdir(initval.plecdir,kymo_output_dir);

kymo_output_dir=[initval.plecdir '\' kymo_output_dir];
outfile_name=[kymo_output_dir '\log.mat'];
%     save(outfile_name,'initval','DNAinfo','tetherpos','start','finish','startI0','finishI0');
save(outfile_name,'initval','DNAinfo');
kymo_file_name=[kymo_output_dir '\kymoghraphs.mat'];
save(kymo_file_name,'KymoRawNick','KymoBgNick','KymoBgBlNick','KymoRawCoil','KymoBgCoil','KymoBgBlCoil','Kymo_DNAdensity','Kymo_pxdensity');
%     end

figure(31)
hdl31=gcf;
savefig(hdl31,[kymo_output_dir '\DNAattachment.fig']);
dlmwrite('Kymograph_DNA.txt',M,'delimiter',' ');
%
% title('Press enter to proceed')
% pause()
% close(31);
%run('MG_pt_detec.m');


