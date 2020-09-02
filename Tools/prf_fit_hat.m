function [prf_hat,prf_residu,fract_hat,fract_res]=prf_fit_hat(prf_kymo,hat,fluo); 
%LP=length(kymo);
corfactor=1;;
prf_hat=hat; 
prf_residu=prf_kymo-corfactor*prf_hat;
sel=find(prf_residu<0);
rel_neg_sum=sum(prf_residu(sel))/sum(hat);
maxiter=10; iter=0;
while (rel_neg_sum<-0.02)&&(iter<maxiter);
iter=iter+1;
prf_residu=prf_kymo-corfactor*prf_hat;
sel=find(prf_residu<0); 
if ~isempty(sel)
    rel_neg_sum=sum(prf_residu(sel))/sum(hat);
    corfactor=corfactor+rel_neg_sum;
    dum=1;
else
    corfactor=corfactor*0.9;
end
end

prf_hat=corfactor*prf_hat;
prf_residu=prf_kymo-prf_hat;

fract_hat=sum(prf_hat)/sum(prf_kymo);
fract_res=1-fract_hat;
if 0
    close all;
    plot(prf_kymo); hold on;
    plot(corfactor*prf_hat);
    plot(prf_residu);
    %[~]=ginput(1); 
    pause(0.1);
end

dum=1;