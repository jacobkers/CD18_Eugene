dircrop='/Users/fabaiwu/Documents/work/drafts/ChrCompactionSegregation/Figure2/material/diffusivity/20161017shorterm/allcropped';

cd(dircrop);
datafiles=dir('*xyco.mat');
for i=1:length(datafiles);
    load(datafiles(i).name);
    plotname=[datafiles(i).name(1:22) 'plot.tif'];

    chrmsdx=0;chrmsdy=0;orimsdx=0;orimsdy=0;termsdx=0;termsdy=0;
%    chrdifx=[];chrdify=[];oridifx=[];oridify=[];terdifx=[];terdify=[];
time=0;
fr=cat(1,1,(1:1:60)');
    for ts=2:61
        time=cat(1,time,(ts-1)*10);
        chrmsdx=cat(1,chrmsdx,chrmsdx(ts-1)+(chrx1(ts)-chrx1(1))^2);
        chrmsdy=cat(1,chrmsdy,chrmsdy(ts-1)+(chry1(ts)-chry1(1))^2);
        orimsdx=cat(1,orimsdx,orimsdx(ts-1)+(orix1(ts)-orix1(1))^2);
        orimsdy=cat(1,orimsdy,orimsdy(ts-1)+(oriy1(ts)-oriy1(1))^2);
        termsdx=cat(1,termsdx,termsdx(ts-1)+(terx1(ts)-terx1(1))^2);
        termsdy=cat(1,termsdy,termsdy(ts-1)+(tery1(ts)-tery1(1))^2);
    end
    chrmsdx=chrmsdx./fr;
    chrdifx=chrmsdx(2:61)./time(2:61)/2;
    chrmsdy=chrmsdy./fr;
    chrdify=chrmsdy(2:61)./time(2:61)/2;
     
    orimsdx=orimsdx./fr;
    oridifx=orimsdx(2:61)./time(2:61)/2;
    orimsdy=orimsdy./fr;
    oridify=orimsdy(2:61)./time(2:61)/2;
    
    termsdx=termsdx./fr;
    terdifx=termsdx(2:61)./time(2:61)/2;
    termsdy=termsdy./fr;
    terdify=termsdy(2:61)./time(2:61)/2;
%% plot data    
    figure(1)
subplot(2,2,1)
    plot(terx1./16,tery1./16,'Color','c')
    hold on;
    plot(orix1./16,oriy1./16,'Color','r')
    plot(chrx1./16,chry1./16,'Color','k')
    xlim([-celll/32 celll/32]);
    ylim([-cellw/32 cellw/32]);
subplot(2,2,2)
    plot(time,termsdx./256,'Color','c');
    hold on;
    plot(time,orimsdx./256,'Color','r');
    plot(time,chrmsdx./256,'Color','k');
    
subplot(2,2,3)
    plot(time(2:61),terdifx,'Color','c');
    hold on;
    plot(time(2:61),oridifx,'Color','r');
    plot(time(2:61),chrdifx,'Color','k');

saveas(gcf,plotname,'tif')
close(gcf);
save(datafiles(i).name,'time','fr','chrmsdx','chrmsdy','chrdifx','chrdify',...
    'orimsdx','orimsdy','oridifx','oridify','termsdx','termsdy','terdifx','terdify','-append');
end
