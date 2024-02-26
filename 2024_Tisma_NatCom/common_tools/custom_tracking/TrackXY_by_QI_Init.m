function  pretracksettings=TrackXY_by_QI_Init(firstim);
            pretracksettings.radialoversampling=2;
            pretracksettings.angularoversampling=0.7;
            pretracksettings.angularoversampling=1;
            pretracksettings.minradius=0;
            pretracksettings.maxradius=25;
            pretracksettings.iterations=10;
            pretracksettings=Build_QI_SamplinggridGrid(pretracksettings) ;
            
function QI=Build_QI_SamplinggridGrid(QI) 
    %Build a radial sampling grid 
    spokesnoperquad=ceil(2*pi*QI.maxradius*QI.angularoversampling/4);
    radbinsno=(QI.maxradius-QI.minradius)*QI.radialoversampling;
    angles=linspace(-pi/4,2*pi-pi/4,spokesnoperquad*4+1)'; 
    angularstep=pi/2/spokesnoperquad;
    angles=angles(1:end-1)+angularstep/2; %to center angles per quadrant
    radbins=linspace(QI.minradius,QI.maxradius,radbinsno);
    [argsgrid,radiigrid]=meshgrid(angles,radbins);
    QI.Y0samplinggrid=(radiigrid.*sin(argsgrid))';
    QI.X0samplinggrid=(radiigrid.*cos(argsgrid))';
    QI.angles=angles;
    QI.radii=radbins;
    