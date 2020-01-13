function info_1=spots0_cleanup_spots(info_1);
%JWJK_C:----[add ABCorC*----------------------------------------------------
%Title: clean spot data from outliers
%Summary: remove spots that are outside a specified range
%References: CD lab, project Eugene Kim, written by Jacob Kers, 2019
%:JWJK_C-----[add ABCorC*---------------------------------------------------
%this function screens and labels  detected  spots
%         info
                    % pos_frameno
                    % pos_X_pix
                    % pos_X_subpix
                    % content_peakvals
                    % content_perspot_est
                    % content_perspot_meas
         info_1.label_OKspot=0*info_1.pos_frameno;         
         Ic=info_1.content_perspot_meas;
  
        %first, get some measures: this one gets the presumed
        %single-condensin counts
        [flag,cleandata]=prf_outlier_flag(Ic,100,0.7,'positive',1); 
        sel=find(flag==1);
        info_1.label_OKspot(sel)=1;
        