function F002a_Process_CellEdgeCheck(initval,actions);
%Track the center for the averaged image and the individual ones  
  
    load(initval.PathResultName);
    
    
    %1 find focus image bf
    [rr,cc,dd]=size(bf_Stack);
    
    
    %2 find focus image rfp
    
    
    drift.x0=zeros(dd,1);
    drift.y0=zeros(dd,1);
    QI=TrackXY_by_QI_Init(bf_pic);
    [x0,y0,~]=TrackXY_by_QI(bf_pic,QI,0);
    QI=CheckEllipticity_by_QI(bf_pic,QI,1);
    for jj=1:dd
        dd-jj+1
        pic=bf_Stack(:,:,jj);   
        [xnw,ynw,~]=TrackXY_by_QI(pic,QI,0); 
        drift.x0(jj)=xnw;
        drift.y0(jj)=ynw;
        pause(0.5);
    end
    %position correction relative to averaged image
    drift.xc=drift.x0-x0;
    drift.yc=drift.y0-y0;
    save(initval.PathResultName,'QI','drift','/append');
    pcolor(bf_pic); colormap bone; shading flat; hold on;

