function trackmap=kym_simulate_twoloops(sho);
%JWJK_C:-------------------------------------------------------------------
%Title: Simulate a kymograph
%Summary: build a simulated, noisy kymograph resembling that of a tether
%containing two loops close to each other. Used to optimize multiple
%tracking.
%Input: none
%Output: kymograph
%References: project Eugene Kim, Jacob Kers 2019
%:JWJK_C-------------------------------------------------------------------
if nargin<1, sho_save=1; end

LT=1000; LX=50; noize=0.2; sigm=1.5;
kymo_a=zeros(LT,LX);  %map        
Taxis=1:LT;
xx_L=round((LX/5)*Taxis/LT+LX/3+1*rand(1,LT));
xx_R=round(xx_L+10+1*randn(1,LT));
peakax=1:LX;       
for ii=1:LT            
   kymo_a(ii,:)=noize*rand(1,LX)+...
           prf_one_gauss_peak(peakax,xx_L(ii),sigm,1)+...
           prf_one_gauss_peak(peakax,xx_R(ii),sigm,1);  
end
trackmap=kymo_a; 
if sho_save==1
    close all
    subplot(1,3,1);
    pcolor(trackmap); shading flat; colormap hot; 
    dlmwrite([pwd, '\Kymograph_DNA.txt'],trackmap);
end
        
        
        