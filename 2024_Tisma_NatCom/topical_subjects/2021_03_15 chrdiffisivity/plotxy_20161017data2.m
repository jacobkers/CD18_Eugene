% This code plots MSD and diffusion rate for single cells

dircrop='/Users/fabaiwu/Documents/work/drafts/ChrCompactionSegregation/Figure2/material/diffusivity/20161017shorterm/allcropped';

cd(dircrop);
datafiles=dir('*xyco.mat');
termsdall=[];orimsdall=[];chrmsdall=[];terdifall=[];oridifall=[];chrdifall=[];celllength=[];
for i=1:length(datafiles)
    load(datafiles(i).name);
    plotname=[datafiles(i).name(1:22) 'plot.tif'];

    chrmsdx=0;chrmsdy=0;orimsdx=0;orimsdy=0;termsdx=0;termsdy=0;
%    chrdifx=[];chrdify=[];oridifx=[];oridify=[];terdifx=[];terdify=[];
time=0;
fr=cat(1,1,(1:1:60)');
    for ts=2:61
        time=cat(1,time,(ts-1)*10);
        chrmsdx=cat(1,chrmsdx,((chrx1(ts)-chrx1(1))^2)/256);
        chrmsdy=cat(1,chrmsdy,((chry1(ts)-chry1(1))^2)/256);
        orimsdx=cat(1,orimsdx,((orix1(ts)-orix1(1))^2)/256);
        orimsdy=cat(1,orimsdy,((oriy1(ts)-oriy1(1))^2)/256);
        termsdx=cat(1,termsdx,((terx1(ts)-terx1(1))^2)/256);
        termsdy=cat(1,termsdy,((tery1(ts)-tery1(1))^2)/256);
    end
    chrdifx=chrmsdx(2:61)./time(2:61)/2;
    chrdify=chrmsdy(2:61)./time(2:61)/2;
    
    oridifx=orimsdx(2:61)./time(2:61)/2;
    oridify=orimsdy(2:61)./time(2:61)/2;
    
    terdifx=termsdx(2:61)./time(2:61)/2;
    terdify=termsdy(2:61)./time(2:61)/2;
    
    celllength=cat(1,celllength,celll/16);
    
    termsdall=cat(2,termsdall,cat(3,termsdx,termsdy,termsdx+termsdy));
    orimsdall=cat(2,orimsdall,cat(3,orimsdx,orimsdy,orimsdx+orimsdy));
    chrmsdall=cat(2,chrmsdall,cat(3,chrmsdx,chrmsdy,chrmsdx+chrmsdy));
    terdifall=cat(2,terdifall,cat(3,terdifx,terdify,(terdifx+terdify)./2));
    oridifall=cat(2,oridifall,cat(3,oridifx,oridify,(oridifx+oridify)./2));
    chrdifall=cat(2,chrdifall,cat(3,chrdifx,chrdify,(chrdifx+chrdify)./2));
    
%% plot data    
%     figure(1)
% subplot(2,2,1)
%     plot(terx1./16,tery1./16,'Color','c')
%     hold on;
%     plot(orix1./16,oriy1./16,'Color','r')
%     plot(chrx1./16,chry1./16,'Color','k')
%     xlim([-celll/32 celll/32]);
%     ylim([-cellw/32 cellw/32]);
%     xlabel('microns'); ylabel('microns');
% subplot(2,2,2)
%     plot(time,termsdx+termsdy,'Color','c');
%     hold on;
%     plot(time,orimsdx+orimsdy,'Color','r');
%     plot(time,chrmsdx+chrmsdy,'Color','k');
%     xlabel('time(sec)'); ylabel('MSD(micron^2)')
% subplot(2,2,3)
%     plot(time(2:61),terdifx+terdify,'Color','c');
%     hold on;
%     plot(time(2:61),oridifx+oridify,'Color','r');
%     plot(time(2:61),chrdifx+chrdify,'Color','k');
%     xlabel('time(sec)');ylabel('D=(MSD/4t)');
% saveas(gcf,plotname,'tif')
% close(gcf);
save(datafiles(i).name,'time','fr','chrmsdx','chrmsdy','chrdifx','chrdify',...
    'orimsdx','orimsdy','oridifx','oridify','termsdx','termsdy','terdifx',...
    'terdify','-append');
end
save('allcells_xy.mat','time','celllength','termsdall','terdifall','orimsdall','oridifall','chrmsdall','chrdifall');
