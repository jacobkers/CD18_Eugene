function A028_Clickit
%collect positions

expname='Figure2Pannel'; 

initval=A001_Initialize_Kymo(expname);
%initval.roilist=[3 26 28]; 

%% if you want to override, redo single ROIs etc.
%initval.roilist=[38]; 

%% get the ROIs to analyze
L_roi=length(initval.roilist);
for ii=1:L_roi

    roistartstop.roino=initval.roilist(ii);  %3 26 %38
    close all;
    Exp=strcat('ROI',num2str(roistartstop.roino));
    SaveName=char(strcat(initval.expi_outpath, Exp,'_clickinfo'));
    datainpath=strcat(initval.expi_inpath,'M', num2str(roistartstop.roino),'\kymo_ImageJ\');       
    source=[ datainpath, initval.kymofile];
    trackmap=dlmread(source);
    [ff,cc]=size(trackmap);
    pcolor(trackmap); colormap(hot); shading flat; hold on;

    %% click the start and stop positions of the loop
    stopit=0; cnt=0;
    while ~stopit
        cnt=cnt+1;
        title('click 3x: 1) pre, 2) start 3) stop  rightclick 3x to end');
        colormap(hot);
        [xx,tt,but]=ginput(3); 
         if but(3)==3, 
                stopit=1;
         else
            roistartstop.pre_t(cnt)=round(tt(1));
            roistartstop.start_x(cnt)=xx(2);
            roistartstop.start_t(cnt)=round(tt(2));
            roistartstop.stop_t(cnt)=round(tt(3));

            %limit time margins
            if roistartstop.pre_t(cnt)<1,roistartstop.pre_t(cnt)=1;end
            if roistartstop.stop_t(cnt)>ff,roistartstop.stop_t(cnt)=ff;end

            %get width of trackbox
            colormap bone
            title('click 2x: L-Rboxwidth');
            [xx,~,~]=ginput(2); 
            roistartstop.trackhalfwidth(cnt)=round((xx(2)-xx(1))/2);
            roistartstop.loopanalysishalfwidth(cnt)=roistartstop.trackhalfwidth(cnt);
         end
    end
    save(SaveName,'roistartstop');
end

