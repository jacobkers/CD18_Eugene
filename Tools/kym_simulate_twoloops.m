function trackmap=kym_simulate_twoloops
%JWJK_C:-------------------------------------------------------------------
%Simulate a kymograph
%Summary: build a simulated, noisy kymograph resembling that of a tether
%containing two loops close to each other. Used to optimize multiple
%tracking.
%Input: none
%Output: kymograph
%References: project Eugene Kim, Jacob Kers 2019
%:JWJK_C-------------------------------------------------------------------

        close all
        initval.BallParkRadius=1;
        initval.PadCurves=1;
        LT=1000;
        LX=50;
        kymo_a=2*rand(LT,LX);  %noisy map
        
        Taxis=1:LT;
        Spotpos1=round((LX/3)*Taxis/LT+LX/3+1*rand(1,LT));
        Spotpos2=round(Spotpos1+10+1*randn(1,LT));
        for ii=1:LT
           for jj=-3:3
           val1=kymo_a(Taxis(ii),Spotpos1(ii)+jj);
           kymo_a(Taxis(ii),Spotpos1(ii)+jj)=val1+(3-abs(jj))^2;

           val2=kymo_a(Taxis(ii),Spotpos2(ii)+jj);
           kymo_a(Taxis(ii),Spotpos2(ii)+jj)=val2+(3-abs(jj))^2;
           end
        end
        trackmap=kymo_a; 
        pcolor(trackmap); shading flat; colormap hot
        dlmwrite([pwd, '\Kymograph_DNA.txt'],trackmap);