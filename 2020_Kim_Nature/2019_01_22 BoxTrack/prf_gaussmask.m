function ar=prf_gaussmask(ar,sigma)
%JWJK_C:----[add ABCorC*----------------------------------------------------
%Summary: This function masks an array with a Gaussian window;
%sigma=1 means edge is on one sigma etc. Used for stable 1D-com
%tracking of a peak mildly centered in a box.
%Input: profile, degree of masking
%Output: masked profile
%References: Jacob Kers, 2019
%:JWJK_C-----[add ABCorC*---------------------------------------------------  
lar=length(ar);
x0=lar/2;
sig=sigma*x0;
ax=linspace(-x0,x0,lar);
ar=ar.*exp(-(ax/sig).^2);