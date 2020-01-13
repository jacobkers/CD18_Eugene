function PP= prf_one_gauss_peak(x,ux,s,normalize)
%JWJK_C:-------------------------------------------------------------------
%Simulate a kymograph
%Summary:  %This is the equation for a normalized1D gaussian peak value one
%Input: axis, mu, sigma, normalize on peak (0) or content (1))
%Output: profile
%References: Jacob Kers 2019
%:JWJK_C-------------------------------------------------------------------
if nargin<4,x=-10:10; ux=0,s=3,normalize=0; end

PP =exp (-(x-ux).^2./(2*s.^2));
if normalize==1
    PP=PP/sum(PP);
end

if nargin<4,close all; plot(PP);end