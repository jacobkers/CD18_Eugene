function prf_res_ft=prf_remove_sloped_background(prf);
    %% Demo mode
    if nargin<1
        prf=prf_make_demo_curves('plectoneme'); 
        prf=prf+linspace(0,1,length(prf)); %slope
    end

section=0.25;
%remove bottom slope of tether profile
    Lp=length(prf); axz=1:Lp;       
    [lo_L,ixL]=min(prf(1:ceil(Lp*section)));
    [lo_R,ixR]=min(prf(end-ceil(Lp*section):end)); ixR=ixR+ceil(Lp-Lp*section+1)-1;
    slopefit=(axz-ixL)*(lo_R-lo_L)/(ixR-ixL)+lo_L;
    %slopefit=polyval(polyfit([ixL ixR],[lo_L lo_R],1),axz);   
    
    prf_res_ft=prf-slopefit;
    
    %% DEMO mode
    if nargin<1      
        plot(prf,'b-'); hold on;
        plot(slopefit,'b-');
        plot(prf_res_ft,'r');       
        legend( 'ori','minima fit',...
                'corrected');
        xlabel('positition, pixel units');
        ylabel('fluorescence intensity, a.u.');
        pause(0.5);        
        [~]=ginput(1);
        hold off;
    end         