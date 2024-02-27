% This is a preliminary code
dircrop='/Users/fabaiwu/Documents/work/drafts/ChrCompactionSegregation/Figure2/material/diffusivity/20161017shorterm/allcropped';

cd(dircrop);
datafiles=dir('*data.mat');
msd=[];
for i=1:length(datafiles);
    load(datafiles(i).name);
    ts=size(celldata,1);
    if ts==61;
        chr1x=chrdata(1,3)-celldata(1,3);
        chr1y=chrdata(1,4)-celldata(1,4);
        chrdif=(chrdata(:,3)-celldata(:,3)-chr1x).^2+(chrdata(:,4)-celldata(:,4)-chr1y).^2;
        focus1x=focidata(1,3,1)-celldata(1,3);
        focus1y=focidata(1,4,1)-celldata(1,4);
        focusdif1=(focidata(:,3,1)-celldata(:,3)-focus1x).^2+(focidata(:,4,1)-celldata(:,4)-focus1y).^2;
        focus2x=focidata(1,3,2)-celldata(1,3);
        focus2y=focidata(1,4,2)-celldata(1,4);
        focusdif2=(focidata(:,3,2)-celldata(:,3)-focus2x).^2+(focidata(:,4,2)-celldata(:,4)-focus2y).^2;
        msd1=[mean(celldata(:,9)) mean(chrdif) mean(focusdif1) mean(focusdif2)];
        msd=cat(1,msd,msd1);
    end
end
numA=msd(:,1)./16;
numB=sqrt(msd(:,2)./256);
numC=sqrt(msd(:,3)./256);
numD=sqrt(msd(:,4)./256);

figure(1);
plot(numA,numD,'LineStyle','none','Marker','d','MarkerSize',6,'Color','b');
xlim([0 15]); ylim([0 2]);
hold on;
%figure(2);
plot(numA,numC,'LineStyle','none','Marker','o','MarkerSize',6,'Color','r');
xlim([0 15]); ylim([0 2]);
%figure(3);
plot(numA,numB,'LineStyle','none','Marker','s','MarkerSize',6,'Color','k');
xlim([0 15]); ylim([0 2]);
