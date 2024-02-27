function pair_table=Get_labelpair_distances(Rfp,Cell,options);
    %JWJK_B:-------------------------------------------------------------------
    %Title:   Classify replication
    %Summary: A single cell image is judged by various distance parameters  
    %Approach: each ori spot is associated with its nearest-neighbour in
    %the same frame. Their common COM point is determined. this point is
    %related to the center-of mass of the cell, and/or the ter position
    %Input: 
    %Cfp:   Rfp with fields: spotY, spotX,spotContent,
    %Cell:  Area:Centroid MajorAxisLength MinorAxisLength Orientation BW Edge 
    %Output: 
        %'dist' table with typical pair distances, 
        %'pair' table with coordinates of nearest-neighbour  pairs
        %[ii xo xo2 xp yo yo2 yp xm ym d_oo]
        %In both cases, indices run by those of the ori spots
    %:JWJK_B-------------------------------------------------------------------
    sep_max=15;
       
    %drop the bad spots
    goodspots=find(Rfp.spotOk);
    Rfp.spotX=Rfp.spotX(goodspots);
    Rfp.spotY=Rfp.spotY(goodspots);
    Rfp.spotContent=Rfp.spotContent(goodspots);
    Rfp.spotOk=Rfp.spotOk(goodspots);
    
    
    xm=Cell.Centroid(1);
    ym=Cell.Centroid(2);
    xxo=Rfp.spotX; 
    yyo=Rfp.spotY; 

    % 1) number of spots;
    N_ori=length(xxo);
    if N_ori<2
        pair_table=[];
    else 
        goodcount=0;
        for ii=1:N_ori   %collect pairs
            cc=0;
            %ori to nearest ori (if N>1) & their shared 'pair center' 
            %are the two  spots also mutually closest?    
            xo=xxo(ii);                  yo=yyo(ii); 
            [xp,yp,d_oo_min, ip]=Find_nearest_neigbour(ii,xxo,yyo);                       
            [~,~,~, ii_nb]=Find_nearest_neigbour(ip,xxo,yyo);            
            mutualpaired=(ii_nb==ii); 
            nearby=d_oo_min<sep_max;
            goodpair=mutualpaired&nearby;
            if goodpair
                cc=cc+1;
                xo2=mean([xo xp]); yo2=mean([yo yp]);
                pair_table(cc)=[cc xo xp xo2 yo yp yo2 xm ym d_oo_min];
            end            
        end
    end 
    
    
    if options.aligncellcenters&~isempty(pair_table);
            pair_table(:,2:4)=pair_table(:,2:4)-xm; pair_table(:,8)=0;
            pair_table(:,5:7)=pair_table(:,5:7)-ym; pair_table(:,9)=0;
    end
    dum=1;
    
    