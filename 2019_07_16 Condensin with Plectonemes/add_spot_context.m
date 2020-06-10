function info_Spot=add_spot_context(info_Spot, info_Otherspot,init)
%this function analyzes spot context: it looks how many spots are nearby a
%spot in the immediate future or the past. 

if nargin<3
 close all;
 pts=3000; frs=1600;
 info_Spot.ori_index=[1];
 info_Spot.pos_roino=[1];
 info_Spot.pos_X_subpix=[50];
 info_Spot.pos_frameno=[frs/2];
 info_Spot.label.OKspot=1;
 info_Otherspot.pos_roino=1*ones(1,pts);  
 info_Otherspot.pos_frameno=frs*rand(1,pts);  % t from 0 to 400
 info_Otherspot.pos_X_subpix=100*rand(1,pts); % x from 0 to 500
 info_Otherspot.label.OKspot=1*ones(1,pts);
end

info_Spot.label.nearto_otherspot_XT=0*info_Spot.label.OKspot;

t_lookaround=200;  %frames to look ahead and back
dx0=5;    %in the next frame, if a spot is within this distance, it is 'connected'
xt_lookahead=0.1;  %this specifies how the x-vicinity is increased per extra frame
%I base it on a quick check: plec-bound cnd wanders some 1 pixel in 10-20 frames

all_roi=info_Spot.pos_roino;
info_Spot.ori_index=1:length(info_Spot.pos_roino);
rois=sort(unique(info_Spot.pos_roino));
counter=0; counts=[];
for ro=1:length(rois)               %for each roi
this_roi=rois(ro);
sel=find(  (info_Spot.pos_roino==this_roi) &...
           (info_Spot.label.OKspot==1));

ori_ix=info_Spot.ori_index(sel);    %the original index
roi_fr=info_Spot.pos_frameno(sel);  %frames
roi_x=info_Spot.pos_X_subpix(sel);  %positions


sel2=find((info_Otherspot.pos_roino==this_roi) &...
           (info_Otherspot.label.OKspot==1));
roi_fr_other=info_Otherspot.pos_frameno(sel2);  %frames of other spots
roi_x_other=info_Otherspot.pos_X_subpix(sel2);  %positions of other spots


Ls=length(roi_fr);
for ii=1:Ls
    counter=counter+1;
    this_ori_idx=ori_ix(ii);
    this_fr=roi_fr(ii);
    this_x= roi_x(ii);
    
    delta_t=abs(roi_fr_other-this_fr);
    delta_x=abs(roi_x_other-this_x);
    nearby=find(...
        (delta_t<=t_lookaround)&...                             %within time limit
        (delta_x<=dx0+delta_t*xt_lookahead));  %within xt-cone, minus self
    if ~isempty(nearby), %add the counts nearby
        neighbour_count=length(nearby);
        mother_index=find(info_Spot.ori_index==this_ori_idx);
        info_Spot.neighbour_count(mother_index)=neighbour_count;
        if neighbour_count>t_lookaround/2
            info_Spot.label.nearto_otherspot_XT(mother_index)=1;
        end
        
        
        %demo plot
        if 0 %ii==round(Ls/2); %nargin <3  %demo mode
            subplot(1,2,1);
            plot(roi_x_other,roi_fr_other,'yo', 'MarkerSize', 3); hold on;
            plot(roi_x_other(nearby),roi_fr_other(nearby),'bo','MarkerSize', 3); hold on;
            plot(this_x,this_fr,'ro', 'MarkerFaceColor', 'r'), hold on;
            legend('all spots','nearby spots','spot of interest');
            xlabel('position');
            ylabel('time');
            [~]=ginput(1); hold off;
        end
    end
end


if 0% nargin<3
      subplot(1,2,2);
      hx=linspace(1,t_lookaround,50);
      data_to_count=info_Spot.neighbour_count;
      hist_neighbour=hist(data_to_count,hx);
      %hist_neighbour(end)=-5;
    
    bar(hx,hist_neighbour,'r');
    ylabel('counts');
    xlabel('#neigbours-t');
    axis tight
    xlim([0 max(hx);]);
     [~]=ginput(1); hold off; close(gcf);
end
end