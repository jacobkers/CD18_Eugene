function kymo_1=Convert_DNA_Kymo(kymo_DNA,levels_DNA);
% levels_DNA
%     level_dark: 1.2870
%     noise_dark: 0.3811
%     tetherlevel: 6.9567
%     level_darktreshold: 2.0493
%     level_tethertreshold: 6.1945
%     level_looptreshold: 7.7190

kymo_1=kymo_DNA-levels_DNA.level_darktreshold;
kymo_1(kymo_1<0)=0;
[rr,cc]=size(kymo_1);
for ii=1:rr
    kymo_1(ii,:)=kymo_1(ii,:)/sum(kymo_1(ii,:))*100;
end
             
if 0
    subplot(1,2,1); 
        pcolor(kymo_DNA); shading flat; colormap bone; hold on;
        title('DNA');
    subplot(1,2,2); 
        pcolor(kymo_1); shading flat; colormap bone; hold on;
        title('perc');
    dum=1;
end