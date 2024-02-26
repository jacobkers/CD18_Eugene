 function [xxip,yyip]=xy_backbone_to_samplinggrid(prf_x,prf_y,shft,alpha);
     %This function builds a sampling grid. it starts with one line XY.
     %That line is duplicated along a specified direction, typically perpendicular to it 
    profilelength=length(prf_x);
    halfwidth=shft;
    tng=tan(alpha/180*pi);
        %expand into sampling grid
    stepvector=shft*linspace(-1,1,2*halfwidth+1);
    xxip=zeros(2*halfwidth+1,profilelength);
    yyip=zeros(2*halfwidth+1,profilelength);
 
    %Then, build a grid from this end use this grid for interpolation; 
     for jj=1:profilelength
            xxip(:,jj)=prf_x(jj)+sin(tng+pi/2)*(stepvector);
            yyip(:,jj)=prf_y(jj)-cos(tng+pi/2)*(stepvector);
     end
%      close all;
%      
%      plot(xxip,yyip,'bo-'); hold on;
%      plot(prf_x,prf_y,'ro-');
%    dum=1;