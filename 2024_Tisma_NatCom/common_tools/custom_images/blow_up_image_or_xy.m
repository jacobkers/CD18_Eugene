function data_out=blow_up_image_or_xy(data_in, blowup, is_BW)
%this function blows up an image, or alternatively, a coordinate column array
%supossed to indicate pixel coordinates (such as contour_x)
%Jacob 2024

if nargin<3, is_BW=0; end
[rr,cc]=size(data_in);
if cc==1 | rr ==1 %1D column data
    is_image=0;
    data_out=data_in*blowup;
else 
    is_image=1;
    mx=max(data_in(:));;
    if blowup>1
        [rr,cc]=size(data_in);
        [XX,YY]=meshgrid(1:cc,1:rr);    
        colax=linspace(1,cc,blowup*cc);
        rowax=linspace(1,rr,blowup*rr);
        [XXip, YYip]=meshgrid(colax,rowax);
        data_in_blw=interp2(XX,YY,data_in,XXip, YYip);
        data_out=JKD2_IM_smoothJK(data_in_blw,blowup);
        back_scaling=mx/max(data_out(:));
        data_out=data_out*back_scaling;
        if is_BW  %back to binary            
            data_out(data_out<0.5)=0;
            data_out(data_out>=0.5)=1;
        end
    else
        data_out=data_in;
    end
end

if 0 && is_image
    close all;
    pcolor(data_out); shading flat, colormap hot; axis equal;
    [~]=ginput(1);
    close(gcf);
end