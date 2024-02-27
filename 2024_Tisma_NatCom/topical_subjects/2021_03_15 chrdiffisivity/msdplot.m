% This is a preliminary code
dircrop='/Users/fabaiwu/Documents/work/drafts/ChrCompactionSegregation/Figure2/material/diffusivity/20160921bn2179007/singlecells';

cd(dircrop);
datafiles=dir('*data.mat');
msd=[];
for i=1:length(datafiles);
    load(datafiles(i).name);
    ts=size(celldata,1);
    if ts==61;
        chr1=chrdata(1,3)-celldata(1,3);
        focus1=focidata(1,3)-celldata(1,3);
        msd1=[celldata(61,7)-celldata(61,5) mean((chrdata(:,3)-celldata(:,3)-chr1).^2) mean((focidata(:,3)-celldata(:,3)-focus1).^2)];
        msd=cat(1,msd,msd1);
    end
end
figure(1);
plot(msd(:,1),msd(:,2),'LineStyle','none','Marker','.','MarkerSize',8);
figure(2);
plot(msd(:,1),msd(:,3),'LineStyle','none','Marker','.','MarkerSize',8);