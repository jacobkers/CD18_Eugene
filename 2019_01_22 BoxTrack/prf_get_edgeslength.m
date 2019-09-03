function [curvestart,curvestop,curveok]=prf_get_edgeslength(prf_res, type_of_profile);
    %function uses 'shaved off' profile to find start and stop
    tresval=0.4;
    Lp=length(prf_res); axz=1:Lp;
        
    [lo_L,ixL]=min(prf_res(1:ceil(Lp/2)));
    [lo_R,ixR]=min(prf_res(ceil(Lp/2):end)); ixR=ixR+ceil(Lp/2)-1;
    slopefit=polyval(polyfit([ixL ixR],[lo_L lo_R],1),axz);
    
    prf_res_ft=prf_res-slopefit;
    
    switch type_of_profile
        case 'tether' %main level, assuming most of profile is tether 
            main_lev=median(prf_res_ft);
            %subplot(2,2,2);
        case 'loop' %%FWHM level, assuming relatively compact loop 
            main_lev=0.5*max(prf_res_ft);
            %subplot(2,2,4);
    end
    lo=min(prf_res_ft);
    sel=find(prf_res_ft>(lo+tresval*(main_lev-lo)));
    
    if ~isempty(sel);
        curvestart=min(sel);     
        curvestop=max(sel);
        curveok=1;
    else
        curvestart=NaN;     
        curvestop=NaN;
        curveok=0;
    end

    if 0      
        plot(prf_res,'b-'); hold on;
        plot(slopefit,'b-');
        plot(prf_res_ft,'r');       
        plot(0*prf_res+main_lev,'r-');
        title(type_of_profile);
        stem(curvestart,prf_res_ft(curvestart),'ro');
        stem(curvestop,prf_res_ft(curvestop),'ro');
        legend( '\fontsize{6} residu','\fontsize{6} minima fit',...
                '\fontsize{6} corrected','\fontsize{6} mean',...
                '\fontsize{6}edges');
        xlabel('positition, pixel units');
        ylabel('fluorescence intensity, a.u.');
        pause(0.5);        
        [~]=ginput(1);
        hold off;
    end         