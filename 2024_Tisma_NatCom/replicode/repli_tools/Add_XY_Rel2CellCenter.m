function LabL=Add_XY_Rel2CellCenter(LabL,Cell);
%This function saves positions relative to the cell center
LabL.spotXc=LabL.spotX-Cell.Centroid(1);
LabL.spotYc=LabL.spotY-Cell.Centroid(2);
dum=1;