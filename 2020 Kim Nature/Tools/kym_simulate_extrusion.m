function [kymo_out_DNA,kymo_out_Cnd]=kym_simulate_extrusion(sho);
%JWJK_C:-------------------------------------------------------------------
%Title: Simulate a kymograph
%Summary: build a simulated, noisy kymograph resembling that of a tether
%extrusion. To be used to  to check our quantification routine.
%Input: none
%Output: kymograph
%References: project Eugene Kim, Jacob Kers 2019
%:JWJK_C-------------------------------------------------------------------
if nargin<1, sho=1; end
LT=1000; LX=60;  sigm=1.5;  
ext_max=30;
Taxis=1:LT;
noize_rel_amplitude_dna=0.03;
noize_rel_amplitude_cnd=0.03;
blank_slate=zeros(LT,LX);  %map 

%% DNA extrusion: 
%starts at 5% at t=200, 
%climbs to 50% over 500 frames, then stays constant at 10%.
%One condensin sits on top
kymo_dna_ext=blank_slate;  %map     
xx_Ext=LX/2+0.2*round((LX/5)*sin(Taxis/LT*6*pi)+1*rand(1,LT));
aa_Ext=(Taxis-200)*ext_max/500; aa_Ext(aa_Ext<0)=0; aa_Ext(aa_Ext>ext_max)=10;
aa_Ext(aa_Ext>0&aa_Ext<5)=5;
peakax=1:LX;    
for ii=1:LT        
   prf=prf_one_gauss_peak(peakax,xx_Ext(ii),sigm,1); 
   prf=prf/sum(prf)*aa_Ext(ii);
   kymo_dna_ext(ii,:)=prf;      
end

%% Extruding condensin
%sits on extrusion location. Content is 20 (arbitrary units)
kymo_cnd_ext=blank_slate;  %map     
for ii=1:LT        
   prf=prf_one_gauss_peak(peakax,xx_Ext(ii),sigm,1);
   isloop=(sum(kymo_dna_ext(ii,:))>0); 
   if isloop
        prf=prf/sum(prf)*10;
        kymo_cnd_ext(ii,:)=prf;   
   end
end

%% Random condensin, not related to loops
%every 50+25 frames, 
%1 condensin lands and leaves at 1/3 and 2/3 ; 
%Content is 15 per condensin (arbitrary units)
kymo_cnd_free=blank_slate;  %map     
for ii=25:50:LT    
        prf=prf_one_gauss_peak(peakax,LX/3,sigm,1)+...;   
            prf_one_gauss_peak(peakax,2*LX/3,sigm,1);
        prf=prf/sum(prf)*2*15;
        kymo_cnd_free(ii,:)=prf;  

end

%% free DNA 'random' plectonemes: 
%every 50 frames, 2 plectonemes, 
%10% content each, %exist for 10 frames, et 1/4h and 3/4th. 
%No condensin involved.
kymo_dna_freeplec=blank_slate;  %map     
for ii=1:50:LT
    for jj=0:9
        prf=prf_one_gauss_peak(peakax,LX/4,sigm,1)+...;   
            prf_one_gauss_peak(peakax,3*LX/4,sigm,1);
        prf=prf/sum(prf)*2*10;
        kymo_dna_freeplec(ii+jj,:)=prf;  
    end
end



%% DNA_tether: 
%a flat background that forms the residu to 100%
kymo_dna_teth=zeros(LT,LX);  %map   
peakposses=(5:0.1:45)/50*LX;
tether_norm=0*peakax;
%make a ''hat'' profile
for pp=1:length(peakposses)
    prf=prf_one_gauss_peak(peakax,peakposses(pp),sigm,1); 
    tether_norm=tether_norm+prf;    
end
tether_norm=tether_norm/sum(tether_norm);
for ii=1:LT        
   used_percentage=sum(kymo_dna_ext(ii,:)+kymo_dna_freeplec(ii,:));
   %used_percentage=sum(kymo_dna_ext(:)+kymo_dna_freeplec(:))/LT;
   rest_perc=100-used_percentage;
   prf=rest_perc*tether_norm;
   kymo_dna_teth(ii,:)=prf;      
end

%% make some noise
noise_cnd=noize_rel_amplitude_cnd*max(kymo_cnd_free(:))*randn(LT,LX);
noise_dna=noize_rel_amplitude_dna*max(kymo_dna_freeplec(:))*randn(LT,LX);


kymo_out_DNA=kymo_dna_ext+kymo_dna_freeplec+kymo_dna_teth+noise_dna; 
kymo_out_Cnd=kymo_cnd_ext+kymo_cnd_free+noise_cnd; 

if sho==1
    close all
    figure;
    subplot(2,1,1);     pcolor(kymo_out_DNA'); shading flat; colormap hot; 
    subplot(2,1,2);     pcolor(kymo_out_Cnd'); shading flat; colormap hot; 
 
    %figure; plot(kymo_out_DNA);
    %dlmwrite([pwd, '\Kymograph_DNA.txt'],kymo_out_DNA);
end
dum=1;     
       
       
        