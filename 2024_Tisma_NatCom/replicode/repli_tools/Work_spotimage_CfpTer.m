function [Cfp, Cfp_pic]=Work_spotimage_CfpTer(Cfp_pic,cellmask,Cell,initval)

%% processing
Cfp_pic=Cfp_pic-median(Cfp_pic(:));    
Cfp_pic=Cfp_pic.*cellmask;
blurPsf=initval.blur_ter;
Cfp_pic=JKD2_IM_smoothJK(Cfp_pic,blurPsf);   
Cfp=Get_MultiSpotProps(Cfp_pic,blurPsf); 
Cfp=Add_XY_Rel2CellCenter(Cfp,Cell);