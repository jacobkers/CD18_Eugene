function [pic,MainPlane]=GetWorkpicFromStack(stack,WorkPicOption);
    stack=double(stack);

    switch WorkPicOption
        case 'MeanProject',
            pic=squeeze(mean(stack,3));
            MainPlane=NaN;
        case 'FocalPlane', 
        %get main plane  (assuming largely planar features)
        stcurve=squeeze(std(std(stack)));
        [~,MainPlane]=max(stcurve);
        pic=double(stack(:,:,MainPlane));  
        case 'First', pic=double(stack(:,:,1));
            MainPlane=1;
    end
    
 