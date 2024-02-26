function [imsh,offset_x]=matrix_pad_verticaledges(im,x);  
JWJK_A:----[add ABCorC*----------------------------------------------------
%Title: padd edges
%Summary: Pads the edge of an image with its mirror as to make it symmetric 
%around a position x
%Input: single-tether roi movies
%Output: spot data of condnesin and plectonemes
%References: CD lab, project Eugene Kim, written by Jacob Kers, 2019
%:JWJK_A-----[add ABCorC*---------------------------------------------------        
[~,cc]=size(im);
offset_x=ceil(x-cc/2);
imsh=im;
if offset_x>0,       
    padstriphi=2*(offset_x);
    padstriplo=1;
    imsh=[im fliplr(im(:,padstriplo:padstriphi))]; 
    imsh=imsh(:,offset_x+1:offset_x+cc);
end
if offset_x<0, 
    padstriplo=cc-(-2*offset_x-1);
    padstriphi=cc;
    imsh=[fliplr(im(:,padstriplo:padstriphi)) im(:,1:padstriphi)]; 
    imsh=imsh(:,-offset_x+1:-offset_x+cc);
end
if 0
    subplot(2,2,1); pcolor(im);
    subplot(2,2,2); pcolor(imsh);
    [~]=ginput(1)
end
            