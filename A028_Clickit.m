function Clickit
%collect positions


expname='Figure2Pannel'; 
%expname='Figure2Pannel'; 

datapath='D:\jkerssemakers\_Data\CD\2018_Eugene\';
exppath=[datapath,'\',expname,'\'];
outpath=strcat(datapath, 'matlabresults\',expname,'\');
textfile='Kymograph_DNA.txt';
roistartstop.roino=38;  %3 26 %38

close all;
Exp=strcat('ROI',num2str(roistartstop.roino));
SaveName=char(strcat(outpath, Exp));
datainpath=strcat(exppath,'M', num2str(roistartstop.roino),'\kymo_ImageJ\');       
source=[ datainpath, textfile];
trackmap=dlmread(source);
[ff,cc]=size(trackmap);

pcolor(trackmap); colormap(hot); shading flat; hold on;

%clickpositions
stopit=0; cnt=0;
while ~stopit
    cnt=cnt+1;
    title('click 2x: start-stop  rightclick 2x to end');
    colormap(hot);
    [xx,tt,but]=ginput(2); 
     if but(2)==3, 
            stopit=1;
     else;
        roistartstop.startx(cnt)=xx(1);
        roistartstop.startt(cnt)=tt(1);
        roistartstop.stopt(cnt)=tt(2);
        if roistartstop.stopt(cnt)>ff,roistartstop.stopt(cnt)=ff;end
        colormap bone
        title('click 2x: L-Rboxwidth');
        [xx,~,~]=ginput(2); 
        roistartstop.trackhalfwidth(cnt)=round((xx(2)-xx(1))/2);
        roistartstop.loopanalysishalfwidth(cnt)=roistartstop.trackhalfwidth(cnt);
     end
end

roistartstop

