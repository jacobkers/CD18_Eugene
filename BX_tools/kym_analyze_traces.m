function looptraces=kym_analyze_traces(looptraces,trackmap,thisroistartstop,ii,initval);  
%JWJK_B:----[add ABCorC*----------------------------------------------------
%Title: Analyze loop traces
%Summary: Position and content measurement per time point
%Approach: For ONE timepoint, all available loops are analyzed. For each
%loop on a time point, it is evaluated if it (still) exists and if so, if
%it can be tracked or if it is in its pre-track low-intensity phase (these
%timepoints are pre-set by user clicks in A028). Tracking is stabilized by
%'lookahead-averaging': a ROI is weighed by the current and immediate
%future timepoints.
%Input: kymograph, a list of loop existnce info and the current time index
%Output: info looptraces
%References: CD lab, project EK, written by JK, 2019
%:JWJK_B-----[add ABCorC*---------------------------------------------------
    %[no xL tL1 tL2 xR tR1 tR2]% 
    Nloops=length(thisroistartstop.start_x);  
   looplist_for_this_frame=[0 0 1];
    
    for iL=1:Nloops %for each loop,perform track analysis in predefined sections of loop
        %% 0) check the loop 'status'
        crp=max([initval.smoothlookahead initval.tracklookahead]);
        loopstandstart=thisroistartstop.pre_t(iL);
        loopwalkstart=thisroistartstop.start_t(iL);
        loopwalkstop=thisroistartstop.stop_t(iL)-crp;
        loopalife=(ii>=loopstandstart&(ii<loopwalkstop));  %does the loop exist on this time?
           
        %% 1) analyze the loop if it exists for this time   
        if loopalife 
            % get an x-position from a kymograph section 'box'
            i_cur=length(looptraces(iL).Lx); %next loop index            
            x_cur=looptraces(iL).Lx(i_cur);
            fr_nxt=looptraces(iL).frame(i_cur)+1;            
            x_nxt=kym_get_next_pos(fr_nxt,iL,ii,x_cur,initval,thisroistartstop,trackmap);
            
            %update positions and times
            looptraces(iL).frame(i_cur+1)=fr_nxt;       
            looptraces(iL).Lx(i_cur+1)=x_nxt;  
            looplist_for_this_frame=[looplist_for_this_frame; [iL, i_cur,x_nxt]];
        end
    end
    %sorting&padding&cleaning
    sel=find(looplist_for_this_frame(:,3)<1);
    if ~isempty(sel)
        looplist_for_this_frame(sel,3)=2;
    end
    [~,idx]=sort(looplist_for_this_frame(:,3));
    looplist_for_this_frame=looplist_for_this_frame(idx,:);
    looplist_for_this_frame=[looplist_for_this_frame; 
                            [0 0 length(trackmap(1,:))]];
    
    
    %% 2, per-loop intensity loop analysis 
    [rr, ~]=size(looplist_for_this_frame);
    N_existingloops=rr-2;
    if N_existingloops>0; %more than just start and stop?
        for jj=1:N_existingloops
            iL=looplist_for_this_frame(jj+1,1);
            xLoop=looplist_for_this_frame(jj+1,3); 
            i_cur=looplist_for_this_frame(jj+1,2); %local index for this loop
            if i_cur==0
                dum=1;
            end
            fL=looptraces(iL).frame(i_cur); 
             
            %% get position, intensity and other basic properties of a newly tracked box; more smoothened
            boxprops.boxhalfwidth=thisroistartstop.loopanalysishalfwidth(iL);
            boxprops.cutlevel=initval.tetherlevel;
            boxprops.lookahead=initval.smoothlookahead;           
            [~,~,prf_box_sm,boxlo,boxhi,box_prf_full]=prf_get_box_intensity(trackmap,xLoop,fL,boxprops);                      
            
            %get edge positions of the loop
            [strt,stp,~]=prf_get_edgeslength(prf_box_sm,'loop');
            loopedge_lft=boxlo-1+strt;
            loopedge_rgt=boxlo-1+stp;
            looptraces(iL).curvestart(i_cur+1)=loopedge_lft;
            looptraces(iL).curvestop(i_cur+1)=loopedge_rgt;
         
            % Intensity measurements
            %build a 'hat-only' section
            prf_box_sm(prf_box_sm<initval.tetherlevel)=initval.tetherlevel;  %padding off-tether
            prf_box_sm=prf_box_sm-min(prf_box_sm);
            
            %build a shaved-off 'residu' tether profile; 
            prf_box_res=box_prf_full;
           
            prf_box_res(boxlo:boxhi)=prf_box_res(boxlo:boxhi)-prf_box_sm; %subtract the loop
    
            if jj==1;  %define also a residu where all hats will be subtracted
                prf_box_res_all=box_prf_full;
                prf_box_full_all=box_prf_full;
            end
            prf_box_res_all(boxlo:boxhi)=prf_box_res_all(boxlo:boxhi)-prf_box_sm; 
            %subtract the current loop  from common tetherprofile            
            
            looptraces=prf_get_regularloop_intensities(looptraces,iL,i_cur,prf_box_sm,prf_box_res,box_prf_full,xLoop);
            looptraces=prf_get_Zloop_intensities(looptraces,iL,i_cur,prf_box_sm,prf_box_res,box_prf_full,loopedge_lft,loopedge_rgt);       
            dum=1;
        end
%         %% 3 run it one time more, now getting tether intensities between neigboring
        %loops from the residu
         for jj=1:N_existingloops
            iL=looplist_for_this_frame(jj+1,1);
            xLoop=round(looplist_for_this_frame(jj+1,3));
            
            x_Nb_lft=round(looplist_for_this_frame(jj,3));
            x_Nb_rgt=round(looplist_for_this_frame(jj+2,3));
            i_cur=looplist_for_this_frame(jj+1,2); %local index for this loop
 
            looptraces=prf_get_inbetween_regularloop_intensities(looptraces,iL,i_cur,prf_box_res_all,prf_box_full_all,xLoop,x_Nb_lft,x_Nb_rgt);
   
            
            dum=1;
          end
        
        
    end
    
    
function x_nxt=kym_get_next_pos(fr_nxt,iL,ii,x_cur,initval,thisroistartstop,trackmap);
%JWJK_C:----[add ABCorC*----------------------------------------------------
%Title: update the new loop position via a kymograph
%Summary: Get a new position from a kymograph using the current one, 
%depending on the loop 'status' : we include a low-intensity
%pre-track where we do get intnesity, but due to the signal weakness, we do
%not try to track.
%References: CDlab, EK, JK, 2019
%:JWJK_C-----[add ABCorC*---------------------------------------------------   
         
%% check the loop 'status'
crp=max([initval.smoothlookahead initval.tracklookahead]);
loopstandstart=thisroistartstop.pre_t(iL);
loopwalkstart=thisroistartstop.start_t(iL);
loopwalkstop=thisroistartstop.stop_t(iL)-crp;
loopalife=(ii>=loopstandstart&(ii<loopwalkstop));  %does the loop exist on this time?
loopstands=(loopalife&(ii<loopwalkstart)); %pre-track section?
loopmoves=(loopalife&(ii>=loopwalkstart)); %track section?

boxprops.boxhalfwidth=thisroistartstop.trackhalfwidth(iL);
boxprops.cutlevel=initval.tetherlevel;
boxprops.lookahead=initval.tracklookahead;
[~,~,box]=prf_get_box_intensity(trackmap,x_cur,fr_nxt,boxprops); %get next box
box(box<initval.tetherlevel)=initval.tetherlevel;  %padding off-tether
box=box-min(box);
mbox=prf_gaussmask(box,0.5);
[com,comc]=prf_get_1D_center_of_mass(mbox);
if loopstands
    x_nxt=x_cur;   %next position; do not update
end
if loopmoves
    x_nxt=x_cur+comc;   % update next position;
end

function looptraces=prf_get_Zloop_intensities(looptraces,iL,i_cur,prf_box_sm,prf_box_res,prf_box_full,loopedge_lft,loopedge_rgt)
%JWJK_C:----[add ABCorC*----------------------------------------------------
%Title: Z-loop intensities
%Summary: obtain intensities corresponding to the Z-loop interpretation;
%this Z loop intensity is made up from two parts: the stem and the hat.
%The stem is the tether section between the estimated edges of the Z-loop.
%The hat is the excess intnesity over the regualr tether intnesity level.   
%Note: a regular loop does not include this stem because then, we assume
%the loop is connected to the tether via one point location)
%References: CDlab, EK, JK, 2019
%:JWJK_C-----[add ABCorC*---------------------------------------------------   
totaltethercounts=sum(prf_box_full);

if loopedge_lft<loopedge_rgt+1
    totalbox_Zhat=sum(prf_box_sm);
    totalbox_Zstem=sum(prf_box_res(loopedge_lft:loopedge_rgt));
    totalbox_Zleft=sum(prf_box_res(1:loopedge_lft-1));
    totalbox_Zright=sum(prf_box_res(loopedge_rgt+1:end));   
else
    totalbox_Zhat=0;
    totalbox_Zstem=0;
    totalbox_Zleft=0;
    totalbox_Zright=0;   
end
looptraces(iL).Zloop.I_left(i_cur+1)=100*totalbox_Zleft/totaltethercounts;
looptraces(iL).Zloop.I_right(i_cur+1)=100*totalbox_Zright/totaltethercounts;           
%Z-loop intensity:
looptraces(iL).Zloop.I_stem(i_cur+1)=100*totalbox_Zstem/totaltethercounts;
looptraces(iL).Zloop.I_hat(i_cur+1)=100*totalbox_Zhat/totaltethercounts;

looptraces(iL).Zloop.I_mid(i_cur+1)= looptraces(iL).Zloop.I_hat(i_cur+1)+...
               looptraces(iL).Zloop.I_stem(i_cur+1);  

looptraces(iL).Zloop.I_checksum(i_cur+1)=...
looptraces(iL).Zloop.I_left(i_cur+1)+...
looptraces(iL).Zloop.I_mid(i_cur+1)+...
looptraces(iL).Zloop.I_right(i_cur+1);     

%note:sum(box_res)+sum(box_sm)=sum(box_full)        
%note: sum(box_res)=sum(box_left)+sum(box_stem)+sum(box_right)     

function looptraces=prf_get_regularloop_intensities(looptraces,iL,i_cur,prf_box_sm,prf_box_res,prf_box_full,xLoop)
%JWJK_C:----[add ABCorC*----------------------------------------------------
%Title:  regular loop intensities
%Summary: obtain intensities corresponding to regular loop interpretation;
%The hat is the excess intensity over the regualr tether intensity level.
%left and right  refer to the full tether residu, including posible other
%loops.
%References: CDlab, EK, JK, 2019
%:JWJK_C-----[add ABCorC*---------------------------------------------------   
totaltethercounts=sum(prf_box_full);

totalbox_hat=sum(prf_box_sm);
xLoop=round(xLoop);
totalbox_left_full=sum(prf_box_res(1:xLoop-1));
totalbox_right_full=sum(prf_box_res(xLoop:end));   

looptraces(iL).Regloop.I_left_full(i_cur+1)=100*totalbox_left_full/totaltethercounts;
looptraces(iL).Regloop.I_right_full(i_cur+1)=100*totalbox_right_full/totaltethercounts;           
looptraces(iL).Regloop.I_hat(i_cur+1)=100*totalbox_hat/totaltethercounts;

looptraces(iL).Regloop.I_checksum(i_cur+1)=...
looptraces(iL).Regloop.I_left_full(i_cur+1)+...
looptraces(iL).Regloop.I_right_full(i_cur+1);     

function looptraces=prf_get_inbetween_regularloop_intensities(looptraces,iL,i_cur,prf_box_res_all,prf_box_full_all,xLoop,x_Nb_lft,x_Nb_rgt);
%JWJK_C:----[add ABCorC*----------------------------------------------------
%Title:  get regular in-between loop intensities
%Summary: obtain intensities corresponding to regular loop interpretation;
%left and right  refer to the tether residu section until the next neighbour
%References: CDlab, EK, JK, 2019
%:JWJK_C-----[add ABCorC*---------------------------------------------------   
totaltethercounts=sum(prf_box_full_all);
totalbox_left_Nb=sum(prf_box_res_all(x_Nb_lft:xLoop-1));
totalbox_right_Nb=sum(prf_box_res_all(xLoop:x_Nb_rgt));   
looptraces(iL).Regloop.I_left_neighbour(i_cur+1)=100*totalbox_left_Nb/totaltethercounts;
looptraces(iL).Regloop.I_right_neighbour(i_cur+1)=100*totalbox_right_Nb/totaltethercounts;           



