function [Rfp,Rfp_pic]=Work_spotimage_RfpOri(Rfp_pic,cellmask,Cell, initval);
    Rfp_pic=Rfp_pic-median(Rfp_pic(:));
    Rfp_pic=Rfp_pic.*cellmask;    
    Rfp=Get_MultiSpotProps(Rfp_pic,initval.Psf_est);   
    Rfp=Screen_spots(Rfp,Rfp_pic,initval,Cell);