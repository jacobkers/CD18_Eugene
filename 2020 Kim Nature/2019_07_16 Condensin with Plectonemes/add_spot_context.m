function info_Spot1=add_spot_context(info_Spot1, info_Spot2,init)
%this function analyzes spot context: it looks how many spots are nearby a
%spot in the immediate future or the past. 

%settings
dx0=5;             %if a next-frame spot is within this distance, it is 'connected'

xt_lookahead_12=0.1;  %this specifies how the x-vicinity is increased per extra frame
                   %I base it on a quick check: plec-bound cnd wanders 
                   %some 1 pixel in 10-20 frames
blinkgap=200;  %frames to look ahead and back
newborngap=10;
lonergap=5;

%% Demo mode
if nargin<3  
 close all;
 frs=1600;
 pts1=1000;
 info_Spot1.ori_index=[1:pts1];
 info_Spot1.pos_roino=1*ones(1,pts1);
 info_Spot1.pos_frameno=[frs/2+frs/3*rand(1,pts1);];
 info_Spot1.pos_X_subpix=[50+40*rand(1,pts1);];
 
 pts2=3000; 
 info_Spot2.pos_roino=1*ones(1,pts2);  
 info_Spot2.pos_frameno=frs*rand(1,pts2);  % t from 0 to 400
 info_Spot2.pos_X_subpix=100*rand(1,pts2); % x from 0 to 500
end
%%

%new fields
info_Spot1.label.nearto_otherspot_XT=0*info_Spot1.pos_roino;
info_Spot1.label.newborn=0*info_Spot1.pos_roino;
info_Spot1.label.loner=0*info_Spot1.pos_roino;

all_roi=info_Spot1.pos_roino;
info_Spot1.ori_index=1:length(info_Spot1.pos_roino);
rois=sort(unique(info_Spot1.pos_roino));
counter=0; counts=[];


%% work spot collections for each roi
for ro=1:length(rois)               
    this_roi=rois(ro);
    if mod(ro,1)==0, disp(strcat('A40_Analyzing spot context roi:',num2str(this_roi),':Exps to work through:',num2str(length(rois)-ro+1)));end  
    %spots
    sel_s1=find(  (info_Spot1.pos_roino==this_roi));
    ori_ix_s1=info_Spot1.ori_index(sel_s1);    %the original index
    roi_fr_s1=info_Spot1.pos_frameno(sel_s1);  %frames
    roi_x_s1=info_Spot1.pos_X_subpix(sel_s1);  %positions
    
    %other spots
    sel2=find((info_Spot2.pos_roino==this_roi));
    roi_fr_s2=info_Spot2.pos_frameno(sel2);  %frames of other spots
    roi_x_s2=info_Spot2.pos_X_subpix(sel2);  %positions of other spots


    Ls=length(roi_fr_s1);
    
    %work all spots1 per time frame
    for ii=1:Ls
        counter=counter+1;
        oridx_i1=ori_ix_s1(ii);
        fr_i1=roi_fr_s1(ii);
        x_i1= roi_x_s1(ii);
        
        %set time and place differences
        delta_t_12=roi_fr_s2-fr_i1;  	%time difference spot 1 to 2
        delta_x_12=abs(roi_x_s2-x_i1);  %place difference spot 1 to 2      
        delta_t_11=roi_fr_s1-fr_i1;      %time difference spot 1 to 1
        delta_x_11=abs(roi_x_s1-x_i1);   %place difference spot 1 to 1 
        mother_index=find(info_Spot1.ori_index==oridx_i1);
        
        % 1) antiblink: look for sufficient nearby spot2 in past, present, future---------------------        
        nearby12_past_future=find(...
            (abs(delta_t_12)<=blinkgap)&...            %within time limit
            (delta_x_12<=dx0+abs(delta_t_12)*xt_lookahead_12));  %within xt-cone, minus self
        if ~isempty(nearby12_past_future) %add the counts nearby
            neighbour_count=length(nearby12_past_future);
            %keep track of the mother index           
            info_Spot1.neighbour_count(mother_index)=neighbour_count;           
            %if more than half of future&past time points have ''other spot'' presence, 
            %the spot-of-interest is deemed ''associated'. Note that a
            %large lookahed bridges blinking gaps of the secondary spot
            if neighbour_count>blinkgap/2
                info_Spot1.label.nearto_otherspot_XT(mother_index)=1;
            end
        end
        %------------------------------------------------------------------
        %2 loners: look for absence of spot1 (self) in near future or past
        nearby11_past_future=find(...
            (abs(delta_t_11)<=lonergap)&...            %within time limit
            (delta_x_11<=dx0+abs(delta_t_11)*xt_lookahead_12));  %within xt-cone, minus self
        if length(nearby11_past_future)<init.t_smooth %nothing nearby (inc. filtering)
            info_Spot1.label.loner(mother_index)=1;
        end
            
        %------------------------------------------------------------------           
         % 3) newborns: look for nearby spot1 (self) only in past--------------------------------   
         nearby11_past=find(...
            ((delta_t_11>-newborngap)&(delta_t_11<0)&...            %within 10 pts time limit
            (delta_x_11<=dx0+abs(delta_t_11)*xt_lookahead_12)));  %within xt-cone, minus self
        if isempty(nearby11_past) %add the counts nearby
            %keep track of the mother index
            info_Spot1.label.newborn(mother_index)=1;
        end
        %------------------------------------------------------------------
        
        
                      
            %% demo plot
            if nargin <3 %ii==round(Ls/2); %nargin <3  %demo mode
                subplot(1,2,1);
                plot(roi_x_s2,roi_fr_s2,'yo', 'MarkerSize', 3); hold on;
                plot(roi_x_s1,roi_fr_s1,'go', 'MarkerFaceColor', 'g','MarkerSize', 1), hold on;
                 plot(x_i1,fr_i1,'ro', 'MarkerFaceColor', 'r','MarkerSize', 6), hold on;
                if ~isempty(nearby12_past_future)
                    plot(roi_x_s2(nearby12_past_future),roi_fr_s2(nearby12_past_future),'bo','MarkerSize', 3); hold on;
                end
                
               
                legend('all spots2','all spots1','spot of interest','nearby spots2');
                xlabel('position');
                ylabel('time');
                [~]=ginput(1); hold off;
            end
        end
end

info_Spot1.label.evaluated=counter;
info_Spot1.label.newborns=sum(1.0*info_Spot1.label.newborn);
info_Spot1.label.loners=sum(1.0*info_Spot1.label.loner);
info_Spot1.label


    if  nargin<3
      subplot(1,2,2);
      hx=linspace(1,blinkgap,50);
      data_to_count=info_Spot1.neighbour_count;
      hist_neighbour=hist(data_to_count,hx);
      bar(hx,hist_neighbour,'r');
      ylabel('counts');
      xlabel('#neigbours-t');
      axis tight
      xlim([0 max(hx);]);
      [~]=ginput(1); hold off; close(gcf);
    end
end