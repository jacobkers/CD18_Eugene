function fitresult=Get_GaussianCrossection(Xsect);
    %Get local Gauss profile (including pre-estimates)
         hw=round((length(Xsect)-1)/2);
         XsectAx=[-hw:hw]';
         
         %normalize etc.
         Ob0=min(Xsect);
         ON0=sum(Xsect-Ob0);
         Xsect_nrm=(Xsect-Ob0)/ON0;
             est.x0=subpix_step(Xsect_nrm);
             est.psf=hw;
             est.b0=0;
             est.N0=1;
        [fitn,YFit]=MLE_One1D_Gaussian_FreePSF(XsectAx,Xsect_nrm,est, 0);
            fitresult.x0=fitn.x0;
            fitresult.N0=fitn.N0*ON0;
            fitresult.b0=fitn.b0*ON0+Ob0;  
            fitresult.psf=fitn.psf;
            fitresult.Fit=YFit*ON0+Ob0;
        
    