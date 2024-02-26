function CleanHoleProps=F007_CleanHoles(SpotProps,HoleProps,im,Psf);
%This function removes holes if they are 'offside'
CleanHoleProps=[];

[LH,~]=size(HoleProps);
[LS,~]=size(SpotProps);
xh=HoleProps(:,3);
yh=HoleProps(:,4);
Ih=HoleProps(:,5);
xs=SpotProps(:,3);
ys=SpotProps(:,4);
Is=SpotProps(:,5);
points=[xs ys Is];

[xm,ym,theta,ecc]=JKD2_XY_calculate2Dmomentpoints(points,0);

        for ii=1:LH          
            HoleX=xh(ii);
            HoleY=yh(ii);
            r_hs=((xs-HoleX).^2+(ys-HoleY).^2).^0.5;
            [dd,ix]=min(r_hs); %nearest spot
            NearestspotX=xs(ix); 
            NearestspotY=ys(ix); 
            R_Spot=((NearestspotX-xm).^2+(NearestspotY-ym).^2).^0.5;
            R_AllSpots=((xs-xm).^2+(ys-ym).^2).^0.5;
            
            
            R_Hole=((HoleX-xm).^2+(HoleY-ym).^2).^0.5; 
            sel=find(R_AllSpots>R_Hole);  
            FarSpotsFract=length(sel)/LH*100;% perc of spots further out
            
            if (R_Hole<R_Spot)|(FarSpotsFract>20)  %if hole is closer than spot, 'inside' 
                CleanHoleProps=[CleanHoleProps; HoleProps(ii,:)];
            end
            dum=1;
        end
