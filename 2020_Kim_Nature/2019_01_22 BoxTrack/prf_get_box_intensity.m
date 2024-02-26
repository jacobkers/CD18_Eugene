function [Iperc,Iperc_excess,box,lox,hix,box_full]=prf_get_box_intensity(trackmap,xx,fr,boxprops);
%JWJK_B:----[add ABCorC*----------------------------------------------------
%Title: obtain the content of a loop and its surroundings.
%Summary: get a local time-averaged local place loop profile and measure
%its (relative) content.
%Approach: a time-place area 'squbox is cut. along the t-axis, a weighted
%time-mask is put before averaging such that the profile from current time
%counts most. A similar procedur is followed for the complete profile; this
%one serves as a reference for 'full DNA content' of the tether.
%Input: kymograph, loop position, frame number, properties of the
%box-to-cut
%Output: content percentages: Irw is the percentage of the box (relative to
%the full tether), Irw_excess is the effective loop intensity; box and
%box_full are the effective time-averaged profiles of local loop and full
%tether, respectively.
%References: CD lab, project Eugene Kim, code by Jacob Kers 2019
%:JWJK_B-----[add ABCorC*---------------------------------------------------

%% 1 obtain the place and time limits
[tt,cc]=size(trackmap);
fr=round(fr);
hf=boxprops.boxhalfwidth;      
lox=round(max([1 xx-hf])); 
hix=round(min([cc xx+hf]));
lofr=round(max([1 fr])); 
hifr=round(min([tt fr+boxprops.lookahead]));

%% 2) build time-averaged profiles
if boxprops.lookahead>0
    squbox=trackmap(lofr:hifr,lox:hix);
    squbox_full=trackmap(lofr:hifr,:);           
    maskline=fliplr(linspace(0.1,1,hifr-lofr+1));
    maskline=maskline/sum(maskline);

    squmask=repmat(maskline',1,hix-lox+1);
    squmask_full=repmat(maskline',1,cc);

    box=sum(squbox.*squmask);
    box_full=sum(squbox_full.*squmask_full);
else
    box=trackmap(fr,lox:hix);
    box_full=trackmap(fr,1:cc);
end
box_full=box_full-min(box_full);

%to avoid dark background issues
box(box<boxprops.cutlevel)=boxprops.cutlevel;

%get the peak excess intensity (above the 'flat tether'level)
Irw=sum(box);
Irw_excess=Irw-boxprops.cutlevel*length(box);

Itotal=sum(box_full);
Iperc=100*Irw/Itotal;
Iperc_excess=100*Irw_excess/Itotal;


    