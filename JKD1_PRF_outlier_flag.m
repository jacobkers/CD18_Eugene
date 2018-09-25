function [flag,cleandata]=JKD1_PRF_Outlier_Flag(data,tolerance,sigchange,how,sho);
% 'Use this section for a Quicksheet'
 %------------------------------------------------------------------
    %this function is meant to find a representative value for a standard
    %deviation in a heavily skewed distribution (typically, flat data with
    % %peaks). It calculates the standard deviation and average the data;
    % Based on these, outliers are determined and excluded for a new calculation
    % of average and SD; this is repeated until sigma does not change anymore too much
    % This is repeated until the new sigma does not change much anymore
    %output: positions of outliers. Jacob Kers 2013 and before
 % 'End of Quicksheet section'

binz=50;

if nargin<5  %For testing/demo purposes
    close all
    data=JK00_DEMODATA_Peaks;
    tolerance=2;
    sigchange=0.7;
    how='positive';
    sho=1;
    plot(data,'o-');
    binz=20;
end

sigma=1E20;            %at start, use a total-upper-limit 
ratio=0;
ld=length(data);
flag=ones(ld,1);  %at start, all points are selected
cleandata=data;
while ratio<sigchange     %if not too much changes anymore; the higher this number the less outliers are peeled off.
    sigma_old=sigma;
    selc=find(flag==1);
    data(flag==1); 
    ls=length(selc);
    av=nanmedian(data(selc));       %since we expect skewed distribution, we use the median iso the mea     
    sigma=nanstd(data(selc));
    ratio=sigma/sigma_old;
    switch how
        case 'positive',  flag=(data-av)<tolerance*sigma;     %adjust outlier flags
        case 'all',  flag=abs(data-av)<tolerance*sigma;     %adjust outlier flags  
    end
    %plot menu------------------  
    if sho==1
        close all;
        cleandata=data(selc); 
        hx=(min(cleandata):(range(cleandata))/binz:max(cleandata));   %make an axis
        sthst=hist(cleandata,hx);
        bar(hx,sthst);
        title('Histogram');
        dum=ginput(1);
        pause(0.5);     
    end
    %---------------------------- 
end
cleandata=data(selc); 
hx=(min(cleandata):(range(cleandata))/binz:max(cleandata));   %make an axis
sthst=hist(cleandata,hx);





