function [info_DNA,info_Cnd]=kym_get_length_and_densities(info_DNA,info_Cnd,kymo_DNA);
%JWJK_C:-------------------------------------------------------------------
%Title: get condnesin density per tether length
%Summary: %This function analyzes the density associated with 
%DNA plectonemes and condensin
%Approach: 
%References: CD lab, project Eugene Kim, written by Jacob Kers, 2019
%:JWJK_C-------------------------------------------------------------------
   
%% get averaged DNA edges
   [rr,cc]=size(kymo_DNA);
   prf_av=nanmedian(kymo_DNA);
   [curvestart,curvestop,curveok]=prf_find_profile_edges(prf_av, 'tether');
   info_DNA.general.curveok=curveok;
   if curveok    
       info_DNA.general.curvestart=curvestart;
       info_DNA.general.curvestop=curvestop;
       info_DNA.general_tetherlength=curvestop-curvestart;
   else
       info_DNA.general_tetherlength=[];
       info_DNA.general.curvestart=[];
       info_DNA.general.curvestop=[];
   end
   
   %% add some classification for the condensin
   isCnd_awayfromDNAedges=(info_Cnd.pos_X_pix>curvestart+3)&...
                      (info_Cnd.pos_X_pix<curvestop-3);
   isPlec_awayfromDNAedges=(info_DNA.pos_X_pix>curvestart+3)&...
                      (info_DNA.pos_X_pix<curvestop-3);
                  
   info_Cnd.classify.awayfromDNAedges=isCnd_awayfromDNAedges;
   info_DNA.classify.awayfromDNAedges=isPlec_awayfromDNAedges;
   
   
   %% get general numbers  for this kymograph
    info_DNA.general_freetetherlength=info_DNA.general_tetherlength-6;
    info_DNA.general_total_tetherlength=info_DNA.general_tetherlength*rr;
    info_DNA.general_total_freetetherlength=info_DNA.general_freetetherlength*rr;   
    info_Cnd.general_total_freetetherlength=info_DNA.general_total_freetetherlength;
    info_Cnd.general_total_number=length(info_Cnd.pos_X_pix);
    
    
    info_Cnd.general_free_number=length(find(...
             (info_Cnd.pos_X_pix>curvestart+3)&...
             (info_Cnd.pos_X_pix<curvestop-3)));
    
     if ~isempty(info_Cnd.general_total_freetetherlength)    
            info_Cnd.general_freedensity=info_Cnd.general_free_number/...
                             info_Cnd.general_total_freetetherlength;
     else
         info_Cnd.general_freedensity=0;
         info_Cnd.general_total_freetetherlength=0;
     end
    dum=1;