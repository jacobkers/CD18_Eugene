function info_Spot=add_persistence_flag(info_Spot,init)
%this function analyes spot context: it looks how many spots are nearby a
%spot in the immediate future or the past
t_lookaround=1000;  %frames to look ahead and back
dx0=2;    %in the next frame, if a spot is within this distance, it is 'connected'
xt_lookahead=0;  %this specifies how the x-vicinity is increased per extra frame
%I base it on a quick check: plec-bound cnd wanders some 1 pixel in 10-20 frames

all_roi=info_Spot.pos_roino;
rois=sort(unique(all_roi));
counter=0; counts=[];
for ro=1:length(rois)
this_roi=rois(ro);
sel=find(all_roi==this_roi);
ori_ix=info_Spot.ori_index(sel); %the original index
roi_fr=info_Spot.pos_frameno(sel);
roi_x=info_Spot.pos_X_subpix(sel);
Ls=length(roi_fr);
for ii=1:Ls
    counter=counter+1;
    this_ori_idx=ori_ix(ii);
    this_fr=roi_fr(ii);
    this_x= roi_x(ii);
    delta_t=abs(roi_fr-this_fr);
    delta_x=abs(roi_x-this_x);
    nearby=find(...
        (delta_t<=t_lookaround)&...                             %within time limit
        (delta_x<=dx0+delta_t*xt_lookahead))-1;  %within xt-cone, minus self
    if ~isempty(nearby), 
        neighbour_count=length(nearby);
        counts(counter)=neighbour_count;
        mother_index=find(info_Spot.ori_index==this_ori_idx);
        info_Spot.neighbour_count(mother_index)=neighbour_count;
    end
end
end