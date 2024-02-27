function results=Localizer_by_envelope(ImA,ImB, showandwait);
%This function ivestigates co-localization of patterns

if nargin<2  
    close all;
    showandwait=1;
    hp=30; noiz=0.1; shf=3;
     xa=5*shf;ya=5*shf;
     xb=5*shf;yb=5*shf;
     
     width=15; 
    [XX,YY]=meshgrid(-hp:hp,-hp:hp);
    RRA=((XX-xa).^2+(YY-ya).^2).^0.5;
    ImA=exp(-(RRA/width).^2)+noiz*rand(2*hp+1);
    RRB1=((XX-xb).^2+(YY-yb).^2).^0.5;
    RRB2=((XX+xb).^2+(YY+yb).^2).^0.5;
    ImB=exp(-(RRB1/width).^2)+noiz*rand(2*hp+1)+...
        1*exp(-(RRB2/width).^2)+noiz*rand(2*hp+1);
end

%% precook
ImA=ImA-min(ImA(:));  %remove background
ImB=ImB-min(ImB(:));

ImA=ImA/max(ImA(:));  %normalize on max
ImB=ImB/max(ImB(:));

ImA_minus_ImB=ImA-ImB;
ImA_minus_ImB(ImA_minus_ImB<0)=0;
ImA_excess=ImA_minus_ImB;
ImA_encompassed=ImA-ImA_excess;


%collect some quantities
results.perc_encompassed=round(sum(ImA_encompassed(:))/sum(ImA(:))*100);

%% plot the results
if showandwait
    figure(2);
    subplot(2,2,1); pcolor(ImA); colormap hot; shading flat;axis equal, axis off;
    title('Channel A');
    subplot(2,2,2); pcolor(ImB); colormap hot; shading flat;axis equal, axis off;
    title('Channel B');
    subplot(2,2,3); pcolor(ImA_excess); colormap hot; shading flat;axis equal, axis off;
    title('Channel A-excess');
    subplot(2,2,4); pcolor(ImA_encompassed); colormap hot; shading flat;axis equal, axis off;
    title('Channel A-encompassed');
    text(0,-0.1*hp*2,[num2str(results.perc_encompassed),'% encompassed']);
     [~]=ginput(1); 
     close(gcf);
end