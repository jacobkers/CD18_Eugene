function results=Localizer_bytresholding(ImA,ImB, showandwait,ABtresholds);
%This function investigates co-localization of patterns

if isnan(ABtresholds) 
    usepresetslimits=0;
else
    usepresetslimits=1;
end

if nargin<2  
    close all;
    showandwait=1;
    hp=100;  shf=hp/7;
     xa=shf;ya=shf; widtha=hp/2; noiza=0.0;
     xb=-2*shf;yb=-2*shf; widthb=hp/3; noizb=0.0;
    [XX,YY]=meshgrid(-hp:hp,-hp:hp);
    RRA=((XX-xa).^2+(YY-ya).^2).^0.5;
    ImA=exp(-(RRA/widtha).^2)+noiza*rand(2*hp+1);
    RRB=((XX-xb).^2+(YY-yb).^2).^0.5;
    ImB=exp(-(RRB/widthb).^2)+noizb*rand(2*hp+1);;
end
ImA=ImA-min(ImA(:));
ImB=ImB-min(ImB(:));

%define treshold lines of various order
if ~usepresetslimits
    [~,tres_A0]=Find_treshold_MD_V2020(ImA);
    [~,tres_B0]=Find_treshold_MD_V2020(ImB);
else
    tres_A0=ABtresholds(1)*max(ImA(:));  %zeroth order
    tres_B0=ABtresholds(2)*max(ImB(:));
end
sel_redzone=find(ImA>tres_A0);
sel_greenzone=find(ImB>tres_B0);
slope=0.25;
sel_whitezone=find((ImA>(tres_A0+ImB*slope))&...
                   (ImB>(tres_B0+ImA*slope)));
sel_yellowzone=find(((ImA<(tres_A0+ImB*slope))|...
                    (ImB<(tres_B0+ImA*slope)))&...
                    ((ImA>tres_A0)&(ImB>tres_B0)));

                lev1=250;
                lev2=220;
                lev3=250;
layer1=0*ImA; 
layer1(sel_redzone)=lev1;  %channel A high
layer2=0*ImB; 
layer2(sel_greenzone)=lev2;  %channel B high
layer3=0*ImB; 
%layer3(sel_whitezone)=lev3; %channel A and B high


%collect some quantities
results.masked_image_A=layer1/lev1;
results.masked_image_B=layer2/lev2;
results.masked_image_A_AND_B=layer3/lev3;
results.area_A=length(find(layer1>0));
results.intensity_A=sum(ImA(find(layer1>0)));
results.area_B=length(find(layer2>0));
results.intensity_B=sum(ImB(find(layer2>0)));
results.area_AB=length(find(layer3>0));
results.intensityA_inAB=sum(ImA(find(layer3>0)));
results.intensityB_inAB=sum(ImB(find(layer3>0)));

results.area_perc_AB_over_A=results.area_AB/results.area_A*100;
results.area_perc_AB_over_B=results.area_AB/results.area_B*100;
results.intensity_percA_AB_over_A=results.intensityA_inAB/results.intensity_A*100;
results.intensity_percB_AB_over_B=results.intensityB_inAB/results.intensity_B*100;

area_im=repmat(0*ImA',1,1,3);

%colorization
area_im(:,:,1)=(layer1');
area_im(:,:,2)=(layer2');
area_im(:,:,3)=(layer3');
area_im=flipud(permute(area_im,[2 1 3]));
area_im=uint8(area_im-1);

%% plot the results
if showandwait
    figure(2);
    subplot(2,3,1); pcolor(ImA); colormap hot; shading flat;axis equal, axis off;
    title('Channel A');
    subplot(2,3,2); pcolor(ImB); colormap hot; shading flat;axis equal, axis off;
    title('Channel B');
    %subplot(2,2,3); plot(ImA(:),ImB(:),'ko', 'MarkerSize',1);  hold on;
    subplot(2,3,3); imshow((area_im)); shading flat;
    title('Overlay tresholds');
    subplot(2,1,2); plot(ImA(sel_redzone),ImB(sel_redzone),'ro', 'MarkerSize',3,'MarkerEdgeColor', 'r','MarkerFaceColor', 'r'); hold on;
        plot(ImA(sel_greenzone),ImB(sel_greenzone),'go', 'MarkerSize',3,'MarkerEdgeColor', 'g','MarkerFaceColor', 'g'); 
        plot(ImA(sel_yellowzone),ImB(sel_yellowzone),'yo', 'MarkerSize',3,'MarkerEdgeColor', 'y','MarkerFaceColor', 'y');hold on;  
        plot(ImA(sel_whitezone),ImB(sel_whitezone),'wo', 'MarkerSize',3,'MarkerEdgeColor', 'k','MarkerFaceColor', 'w');    
        xlabel('intensity A');
        ylabel('intensity B');
        legend('A high','B high','A&B high', 'A,B and ratios high')
        legend('Location', 'Eastoutside');   
     [~]=ginput(1);
   
    outname=strcat('C:\Users\jkerssemakers\Dropbox\CD_recent/temp.tif');
    imout=uint8(area_im-1);
    imwrite(imout,outname); 
     close(gcf);   
end