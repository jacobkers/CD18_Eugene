function cell_basics_num=Get_BasicCellGeometry(Chromosome);
    %Collect average properties of cells:
    %Cromosome radius; peaks
            
    Chro_Rp_av = nanmean(Chromosome.PolarContourMax);
    Chro_Rp_st = nanstd(Chromosome.PolarContourMax);
    Chro_Rp_mx = nanmax(Chromosome.PolarContourMax);
    Chro_Rp_mn = nanmin(Chromosome.PolarContourMax);
    
    
    %Cromosome radius; outer contour
    Chro_Ro_av = nanmean(Chromosome.PolarContourEdge);
    Chro_Ro_st = nanstd(Chromosome.PolarContourEdge);
    Chro_Ro_mx = nanmax(Chromosome.PolarContourEdge);
    Chro_Ro_mn = nanmin(Chromosome.PolarContourEdge);
    
    %Cromosome radius of gyration
    Chro_Rg = Chromosome.RadGyr;
    
    %Cromosome FWHM
    Chro_W_av = nanmean(Chromosome.PolarContourFWHM);
    Chro_W_st = nanstd(Chromosome.PolarContourFWHM);
    Chro_W_mx = nanmax(Chromosome.PolarContourFWHM);
    Chro_W_mn = nanmin(Chromosome.PolarContourFWHM);
    
    %chromosome length
    Chro_Lp=Chromosome.TotalMaxPeakLength;
    Chro_Lo=Chromosome.TotalContourLength;
        
    % cell outline    
    CellWall_R_av=nanmean(Chromosome.RadialCellEdge);
    CellWall_R_st=nanstd(Chromosome.RadialCellEdge);
    CellWall_R_mx=nanmax(Chromosome.RadialCellEdge);
    CellWall_R_mn=nanmin(Chromosome.RadialCellEdge);
    
    % get total lengths cell wall
    CWX=Chromosome.CartesianCellwallX;
    CWY=Chromosome.CartesianCellwallY;                             
    CellWall_L=nansum(((CWX(2:end)-CWX(1:end-1)).^2+(CWY(2:end)-CWY(1:end-1)).^2).^0.5);    
              
    %numeric output
    cell_basics_num=...
       [Chro_Lp ...
        Chro_Lo ...
        Chro_Rp_av Chro_Rp_st Chro_Rp_mx Chro_Rp_mn ...
        Chro_Ro_av Chro_Ro_st Chro_Ro_mx Chro_Ro_mn Chro_Rg ...
        Chro_W_av Chro_W_st Chro_W_mx Chro_W_mn ...
        CellWall_L CellWall_R_av CellWall_R_st CellWall_R_mx CellWall_R_mn]; 
    