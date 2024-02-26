function F002_Process_BrightField_and_Cytosol(initval,actions);
%Track the center for the averaged image and the individual ones  
  
    load(initval.PathResultName);    
    QI=TrackXY_by_QI_Init(bf_pic);
    [x0,y0,~]=TrackXY_by_QI(bf_pic,QI,1);
    QI=CheckEllipticity_by_QI(bf_pic,QI,1);

    save(initval.PathResultName,'QI','/append');
    pcolor(bf_pic); colormap bone; shading flat; hold on;
    
    
    
    %now, some post-processing
    Xsamplinggrid=QI.X0samplinggrid+QI.xpos;
    Ysamplinggrid=QI.Y0samplinggrid+QI.ypos;
    allprofiles_bf=(interp2(bf_pic,Xsamplinggrid,Ysamplinggrid,'NaN'));
    allprofiles_rfp=(interp2(rfp_pic,Xsamplinggrid,Ysamplinggrid,'NaN'));
    
    
    [aa,radii]=size(allprofiles_bf);
    
    
    %alignment of plots and curve
    SamplingAngles=QI.angles; 
    EdgeAngles=QI.CellEdgesVsangle(:,1);
    startangle=SamplingAngles(1);
    [EdgeAngles,ix]=sort(mod(EdgeAngles+startangle,2*pi));
    MappingIndex=(1:aa)'; 
    MappingIndex=MappingIndex(ix);
    [MappingIndex,ix2]=sort(MappingIndex);
    CellBfEdge=QI.radialoversampling*QI.CellEdgesVsangle(ix2,2)
    
    
    
    
    
    figure(4);
    subplot(1,2,1); 
    pcolor(allprofiles_bf); colormap bone; shading flat; hold on;
    plot(CellBfEdge, MappingIndex, 'r-', 'Linewidth', 2);
    title('radial map of brightfield with bf-edge overlay');
    xlabel('radial sampling steps');
    ylabel('angular sampling steps');
    subplot(1,2,2); 
    pcolor(allprofiles_rfp); colormap bone; shading flat; hold on;
    plot(CellBfEdge, MappingIndex, 'r-', 'Linewidth', 2);
    title('radial map of rfp-cytosol with bf-edge overlay');
    xlabel('radial sampling steps');
    ylabel('angular sampling steps');
    
    
    save(initval.PathResultName, 'allprofiles_bf', ...
                                 'allprofiles_rfp', ...
                                 'CellBfEdge',...
                                 '/append');
    
    dum=1;