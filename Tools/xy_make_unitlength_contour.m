function [X,Y]=xy_make_unitlength_contour(X,Y,reps)
        ye=Y(end); xe=X(end);
        for cc=1:reps
            %measure contour length; determine average distance between points
            CL=round(sum(((X(2:end)-X(1:end-1)).^2+...
                    (Y(2:end)-Y(1:end-1)).^2).^0.5));  %approximate contour length          
            [X,Y]=xy_get_smooth_xyline(X,Y,CL,5);
            Y=[Y'; ye]; X=[X'; xe];
        end