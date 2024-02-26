function idx=Get_GuidedEdge(Profile,startpos)
%This function finds a certain value, prefererentially in a pre-specified region
%JacobKers2016
edgeval=0.3;  %fraction of nearby mask

if nargin<3    
    if nargin<1  %no data
        close all;
        ax=-50:50;
        Profile=exp(-0.5*((ax+40)/15).^2); 
        plot(Profile,'o-');
        [~,startpos]=max(Profile);
    end
    Ld=length(Profile);
    MaskCenter=Ld/2;
    MaskWidth=Ld/4;
end
Profile=Profile-min(Profile);
Profile=Profile/Profile(startpos);  %local normalization
Ld=length(Profile);
ax=(1:Ld);

stopit=0;
ix1=startpos;
difcurve=abs(Profile-edgeval);
while~stopit
     ix2=ix+1;
     slope=difcurve(ix2)-difcurve(ix1);
     if slope>=0, stopit=1
     
     
end