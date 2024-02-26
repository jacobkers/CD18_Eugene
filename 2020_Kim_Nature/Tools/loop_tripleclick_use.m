function kym_tripleclick_use
%JWJK_C:----[add ABCorC*----------------------------------------------------
%Title: click a triple-point loop and save the trace
%Summary: convert from tif to txt
%:JWJK_C------[add ABCorC*---------------------------------------------------
%Initialize section--------------------------------------------------------
%file handling: setup general paths 
close all;
inpth='C:\Users\jkerssemakers\CD_Data_in\2018_Eugene\2019_10_14 slippage\';
outpth='C:\Users\jkerssemakers\Dropbox\CD_Data_out\2018_Eugene\2019_10_14 slippage\';
expi=2;
switch expi
    case 1
        label='tether1';
        source=strcat(inpth,'190831_173555_motor_slippage-1_smz.tif');
        target=strcat(outpth,label, '_results.txt');
        clicks=dlmread([outpth,'tripleclicks_190831_173555_motor_slippage-1.txt']);
    case 2
        label='tether2';
        source=strcat(inpth,'slippage_2nd.tif');
        target=strcat(outpth,label,'_results.txt');
        clicks=dlmread([outpth,'tripleclicks_2nd.txt']);
        
end

%build axis
fr_ax=(clicks(1,1):clicks(end,1));
xx=interp1(clicks(:,1),clicks(:,2),fr_ax,'linear');
yy=interp1(clicks(:,1),clicks(:,3),fr_ax,'linear');

FF=length(fr_ax);



firstplane=double(imread(source,'Index',1));
info = imfinfo(source);  
[rr,cc]=size(firstplane);     

%1)  all images
results=[];
for idx=1:FF
    fr_idx=fr_ax(idx)
    pic=double(imread(source,'Index',fr_idx));
    if 0
        figure (1); pcolor(pic); shading flat; colormap jet; hold on;
            title(num2str(fr_idx))
        plot(xx(idx),yy(idx),'o', 'MarkerFaceColor', 'w'); hold off;
        pause(0.1);
        dum=1;
    end
    
    %% intensity analysis
    if 0
    I_lp=get_intensity(expi,pic,xx(idx),yy(idx),'loop');
    I_ms=get_intensity(expi,pic,xx(idx),yy(idx),'motorside');
    I_as=get_intensity(expi,pic,xx(idx),yy(idx),'anchorside');
    I_aa=I_lp+I_ms+I_as;
    Gen_lp=I_lp/I_aa*100;
    Gen_ms=I_ms/I_aa*100;
    Gen_as=I_as/I_aa*100;
    results=[results; [fr_idx Gen_lp Gen_ms Gen_as]];
    end
    %% geometry analysis
    L_ms=get_length(expi,pic,xx(idx),yy(idx),'motorside');
    L_as=get_length(expi,pic,xx(idx),yy(idx),'anchorside');
    results=[results; [fr_idx L_ms L_as]];
    
end
close all;
plot(0.1*results(:,1),0.125*results(:,2:end), 'Linewidth', 2);
xlabel('time, seconds');
ylabel('distance, microns');
%xlim([100 140]);
%legend('loop','motor side','anchor side')
legend('motor side','anchor side')
dlmwrite(target, results);
saveas(gcf,strcat(outpth,label,'_plot.fig'));
saveas(gcf,strcat(outpth,label,'_plot.jpg'), 'jpeg');
%save(target, results);


function I_sum=get_intensity(expi,pic,x,y,fromwhat);
%get areas from bounding box
switch expi
    case 1  %first loop
            switch fromwhat    
                case 'loop' 
                    exclude=0;
                    boundingbox=[-77 37];  %hix lox hiy loy    
                    lox=round(x+boundingbox(1)); hix=round(x-exclude);
                    loy=round(y+exclude); hiy=round(y+boundingbox(2));
                case 'motorside' 
                    exclude=1;
                    boundingbox=[51 -26];  %hix lox hiy loy    
                    lox=round(x+exclude); hix=round(x+boundingbox(1));
                    loy=round(y+boundingbox(2)); hiy=round(y-exclude);
               case 'anchorside' 
                    exclude=1;
                    boundingbox=[17 42];  %hix lox hiy loy    
                    lox=round(x+exclude); hix=round(x+boundingbox(1));
                    loy=round(y+exclude); hiy=round(y+boundingbox(2));
              case 'all' 
                    boundingbox=[-77 51 -26 42];  %hix lox hiy loy    
                    lox=round(x+boundingbox(1)); hix=round(x+boundingbox(2));
                    loy=round(y+boundingbox(3)); hiy=round(y+boundingbox(4));
            end
    case 2  %second loop
            switch fromwhat    
                case 'loop' 
                    exclude=0;
                    boundingbox=[-22 12];  %hix lox hiy loy    
                    lox=round(x+boundingbox(1)); hix=round(x-exclude);
                    loy=round(y+exclude); hiy=round(y+boundingbox(2));
                case 'motorside' 
                    exclude=2;
                    boundingbox=[14 12];  %hix lox hiy loy    
                    lox=round(x+exclude); hix=round(x+boundingbox(1));
                    loy=round(y+exclude); hiy=round(y+boundingbox(2));
               case 'anchorside' 
                    exclude=2;
                    boundingbox=[-22 -30];  %hix lox hiy loy    
                    lox=round(x+boundingbox(1)); hix=round(x-exclude);
                    loy=round(y+boundingbox(2)); hiy=round(y-exclude);
              case 'all' 
                    boundingbox=[-22 14 -30 12];  %hix lox hiy loy    
                    lox=round(x+boundingbox(1)); hix=round(x+boundingbox(2));
                    loy=round(y+boundingbox(3)); hiy=round(y+boundingbox(4));
            end
end
area_sel=pic(loy:hiy,lox:hix);
area_sel=area_sel-median(area_sel(:));
I_sum=sum(area_sel(:));


function L=get_length(expi,pic,x,y,fromwhat);
%get areas from bounding box
switch expi
    case 1  %first loop
            switch fromwhat    
                case 'motorside', tetherx=140; tethery=30;  
               case 'anchorside', tetherx=124; tethery=78; 
                    
            end
    case 2  %second loop
           switch fromwhat    
                case 'motorside', tetherx=61; tethery=47;  
               case 'anchorside', tetherx=53; tethery=17; 
    end
end
if 0
    figure (2); pcolor(pic); shading flat; colormap jet; hold on;     
        plot(x,y,'o', 'MarkerFaceColor', 'w');
        plot(tetherx,tethery,'o', 'MarkerFaceColor', 'w'); hold off;
        pause(0.1);
end

L=((x-tetherx).^2+(y-tethery).^2).^0.5;


