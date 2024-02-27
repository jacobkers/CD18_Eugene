function chro_pic=Remove_Background(chro_pic,howtodoit);

switch howtodoit
    case 'Median', RealMinVal=nanmedian(chro_pic(chro_pic~=0));
    case 'Min', RealMinVal=nanmin(chro_pic(chro_pic~=0));
end
chro_pic(chro_pic==0)=RealMinVal;
chro_pic=chro_pic-nanmedian(chro_pic(:));
