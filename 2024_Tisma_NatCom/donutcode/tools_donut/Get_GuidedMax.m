function idx=Get_GuidedMax(Profile,MaskCenter,MaskWidth)
    %JWJK_C:
    %This function finds a maximum along a profile, but prefererentially in 
    %a pre-specified region. This is done by iterative masking. 
    %:JWJK_C
if nargin<3    
    if nargin<1  %no data
        ax=-50:50;
        Profile=exp(-0.5*((ax-40)/15).^2); 
        plot(Profile,'o-');
    end
    Ld=length(Profile);
    MaskCenter=Ld/2;
    MaskWidth=Ld/4;
end
Ld=length(Profile);
ax=(1:Ld);

for ii=1:5
    Mask=exp(-0.5*((ax-MaskCenter)/MaskWidth).^2); 
    MaskedProfile=Mask.*Profile;
    [~,idx]=max(MaskedProfile);
    MaskCenter=idx;
end