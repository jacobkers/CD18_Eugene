function [prf_post_hat2,prf_loop,fract_hat,fract_loop]=prf_fit_hat(prf_kymo,prf_pre_hat,fluo,init,tempbreak, method); 
%This function splits a tether plus loops profile in its two contributing
%profiles: a 'hat' profile and a 'loop' profile. There are various options
%to find the best division

switch method
    case 'push'
        %a range of hat-curves is detrmined, ranging around an expected
        %amplitude. Then, atresholding routine finds the point that the
        %amplitude of the hat curve gets too large (i.e., the number of
        %erroneously subtracted values rises strongly)
        max_iter=20;
        corfactors=linspace(0.3, 2,max_iter); 
        stp=median(diff(corfactors));
        hat_matrix=corfactors'*prf_pre_hat;
        prf_matrix=repmat(prf_kymo,max_iter,1);
        dif_matrix=prf_matrix-hat_matrix;
        count_curve=sum(1.0*dif_matrix'<0);  %count only overshoots
        dif_matrix(dif_matrix>0)=0;           %sum only overshoots       
        overshootsum=-nansum((dif_matrix'))./sum(prf_matrix'); 
        %find a treshold where many points get overcorrected
        [~,corfactor,~]=prf_find_treshold_MD2020(corfactors,overshootsum,0);
        corfactor=corfactor(1);   
        prf_post_hat1=corfactor*prf_pre_hat;
        prf_loop=prf_kymo-prf_post_hat1;     
    case 'direct_subtract'  %simplest case: constant subtraction
        prf_post_hat1=prf_pre_hat;        
    case 'medians_align'    %ad
        med1=median(prf_pre_hat);
        med2=median(prf_kymo);
        corfactor=med2/med1
        prf_post_hat1=prf_pre_hat*corfactor;
end
prf_loop=prf_kymo-prf_post_hat1;  %split with best hat estimate
%post cooking. Important while doing this is, the original content of the kymograph is
%preserved (for example, if this was 100%, sum of hat and loops should be just that)
prf_loop(prf_loop<0)=0; %shave off 'pits'
prf_post_hat2=prf_kymo-prf_loop;
%checksums:
sumhat=sum(prf_post_hat2);
sumloop=sum(prf_loop);
sumori=sum(prf_kymo);
checkdif=sumori-sumhat-sumloop; %should be zero

fract_hat=sumhat/sumori;
fract_loop=sumloop/sumori;

if tempbreak
    corfactor
    close all;
    plot(prf_kymo); hold on;   
    plot(prf_post_hat2);
    plot(prf_loop, 'Linewidth', 2); hold on;
    legend('kymo','hat','loop');
    [~]=ginput(1); 
    pause(0.1);
end

dum=1;

% OLD CORFACTOR CODE
% prf_hat=(smooth(hat',ceil(2*init.psf_est)))'; 
%         prf_loop=prf_kymo-corfactor*prf_hat;
%         sel=find(prf_loop<0);
%         overshoot=-sum(prf_loop(sel))/sum(prf_kymo);
%         
%         overshoot_mx=0.1*length(prf_loop);
%         overshoot=length(sel);
%         %if 'hat' gets too large, overshoot gets positive
% 
%         %iterate to push the hat profile, properly scaled, from below against the full
%         %profile. Note the small overshoot allowed. The positive difference is the loops. 
%         
%         step=0.05;
%         iter=0;       
%         while (overshoot<=overshoot_mx)&&(iter<max_iter);
%             iter=iter+1;
%             prf_loop=prf_kymo-corfactor*prf_hat;
%             sel=find(prf_loop<0); %count overshoot pixels
%             overshoot=length(sel);              
%             corfactor=corfactor+step  %increase correction
%             dum=1;
%         end