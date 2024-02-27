function Get_noise_from_raw_image
close all;
clear all
for ii=1:11
    
 switch 2
     case 2
         %xy 1:11
         xy=ii;
        inpth='D:\jkerssemakers\_Data\CD\2016_Sandro\2016_11_23_2179_Wildtype_Zcheck\tiffraw\';
        in_name=['20161123bn2179a22001xy',num2str(xy,'%02.0f'),'c2z09']; 
         out_name=['20161123bn2179a22001','xyall','%02.0f','c2z09'] ;
     case 2
         %noizes=[0  0.02  0.05  0.2 0.5]
        rep=ii;
        inpth='D:\jkerssemakers\_Recent\CD\BN_CD16_Sandro\Matlabcode\2018_09_TestingClusters\Results_gen_2_fullthick\';
        in_name=['Sim_CL06DW0.1Noise0.05Repeat',num2str(rep,'%02.0f'),'psf_struct5.0psf_opt2.5'];
        out_name=['Sim_CL06DW0.1Noise0.05','_allreps','psf_struct5.0psf_opt2.5'];
 end
 load_tif=strcat(inpth,in_name,'.tif'); 
 im=double(imread(load_tif));  
 
 edgedata=[im(1:2,1:end) ;im(end-1:end,1:end)];
 [~,edgedata]=JKD1_PRF_outlier_flag(edgedata(:),3,0.7,'positive',0);
 
 noise(ii)=std(edgedata);
 level(ii)=mean(edgedata);
 
 vals=sort(im(:));
 mn=vals(1); 
 mx=vals(end);
 rng=vals(end)-vals(1);
 dix=find((vals>mx-0.6*rng)&(vals<mx-0.4*rng));
 donutvals_sort=vals(dix);
 LD=length(donutvals_sort);
 axz=(1:LD)';
 
 sel=find(im>(mx-0.6*rng)); 
 prms=polyfit(axz,donutvals_sort,1); %linear fit to extrapolate drop-off
 Lfit=axz*prms(1)+prms(2); 
 backboneval(ii)=Lfit(end);
 
if ii==4
    figure;
     subplot(2,2,1);
        pcolor(im); shading flat;
    subplot(2,2,2);
        plot(edgedata);
        title('edge values');
        xlabel('pixel index');
        ylabel('value, a.u.');
    subplot(2,2,3);
        plot(im(sel));
        title('donut values');
        xlabel('pixel index');
        ylabel('value, a.u.');

     subplot(2,2,4);
     plot(axz,donutvals_sort); hold on;
     plot(axz,Lfit,'r-');
        title('donut ridge fit');
        xlabel('pixel index,valaue sorted');
        ylabel('value, a.u.');
end
 
NSratio(ii)=noise(ii)/backboneval(ii);
end
noise'
av_noise=mean(noise)
backboneval'
av_backboneval=mean(backboneval)
NSratio'
av_NSratio=mean(NSratio)
std_NSratio=std(NSratio)

save([inpth,out_name,'_measured_noise.mat'],'noise', 'av_noise', 'backboneval','av_backboneval','NSratio','av_NSratio', 'std_NSratio');
 
 