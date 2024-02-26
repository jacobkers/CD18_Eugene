function [idx,stp, avl, avr,varL,varR, rank1,rank2]=Splitfast(r,initval);               %
    %this function also adresses a one-dim array 'r' and determines the best step-fit there
    %To save (a lot of) time, functions like 'mean' are avoided and results of the former
    %step-fit are used. 
    steprepulsion=0.1;
     w=length(r);
    sk=2;          %skip edges: steps are supposed to stay more than this away from each other
    
    if w>2*sk+1       %plateaus of size '2sk+1' and smaller are not split
        %1) The first fit starts fully at the left
		Chisq=(1:w)*0;  %was (1:w-1)  
		AvL=sum(r(1:sk))/(sk);                             
        AvR=sum(r(sk+1:w))/(w-sk); 
        AvL_old=AvL;
        AvR_old=AvR;
        sumVarL=sum(r(1:sk).^2)-(sk)*AvL^2;;
        sumVarR=sum(r(sk+1:w).^2)-(w-sk)*AvR^2;
        
        %2) Then, all locations are tried for a step-fit; 
        %a chi-square array is built
        for t=sk+1:w-sk
            %r(t) is the value that is going to be transferred from the
            %left to the right plateau; the new index indicating the step
            %(always prior to change of plateau) will thus be ix=t-1
            ix=t-1;      
            AvL=(AvL_old*ix+r(ix+1))/(ix+1);        %left average is adapted, using one new value
            AvR=(AvR_old*(w-ix)-r(ix+1))/(w-ix-1);  %right average is adapted, using one new value 
            sumVarL=sumVarL+r(ix+1)^2 +(ix)*AvL_old^2-(ix+1)*AvL^2;     %left variance value is adapted from its former value
            sumVarR=sumVarR-r(ix+1)^2 +(w-ix)*AvR_old^2-(w-ix-1)*AvR^2; %right variance value is adapted from its former value           
            asy=steprepulsion*(1/(1-1/(t))+1/(1-1/(w-t+1)));     %edge fit compensator
            
            Chisq(ix+1)=(sumVarL+sumVarR)/w *(1+asy);       %weighted average of two chisquares, compensated for edge effect 
            AvL_old=AvL;                            %old variance is stored, left plateau
            AvR_old=AvR;                            %old variance is stored, right plateau
        end
        [rank1,idx]=min(Chisq(sk+1:w-sk)); %the best fit is that with the lowest Chi-square
        if 0
            idx=round(w/2);
            rank1=Chisq(idx);
        end
        idx=idx+2; 
        if 0
            subplot(2,1,1); plot(r(sk+1:w-sk)); subplot(2,1,2); plot(Chisq(sk+1:w-sk),'r'); 
            [a,b,c]=ginput(1); close(gcf);
        end
          
         %3) determine final fit properties   
         avl=sum(r(1:idx))/(idx);                                         
         avr=sum(r(idx+1:w))/(w-idx);
         varL=(sum((r(1:idx)-avl).^2))/(idx);     
         varR=(sum((r(idx+1:w)-avr).^2))/(w-idx);
         stp=avr-avl;
         rank2=stp^2/(1/(idx-1)+1/(w-idx));  
            %quantity expresing accuracy step (squared)
            %bases on expected relative error in step     
     else
         idx=ceil(w/2);
         stp=0; 
        avl=sum(r)/w; 
        avr=avl;
        varL=sum(r.^2)-w*avl^2; ;
        varR=varL ;
        rank1=0;
        rank2=0;
         %low-limit cases-------------------------------------
         
         switch w 
             case 3                             %plateau with only 3 values
                a=(r(2)+r(3))/2-r(1); 
                b=r(3)-(r(1)+r(2))/2;
                cL=[r(1) , (r(1)+r(2))/2];
                cR=[(r(2)+r(3))/2 , r(3)]; 
                [stp,idx]=max([a b]);
                avl=cL(idx);            avr=cR(idx);;
            case 2                              %plateau with only 2 values
                idx=1; stp=r(2)-r(1); avl=r(1); avr=r(2); 
            case 1
               % idx=1; stp=0; avl=r(1); avr=r(1);
        end
    
    end