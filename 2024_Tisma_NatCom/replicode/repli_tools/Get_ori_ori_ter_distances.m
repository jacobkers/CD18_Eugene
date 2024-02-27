function ori_ori_ter_table=Get_ori_ori_ter_distances(Cfp,Rfp,Cell);
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
        %[ii xo xo2 xp yo yo2 yp xm ym]
        %In both cases, indices run by those of the ori spots
    %:JWJK_B-------------------------------------------------------------------
    xm=Cell.Centroid(1);
    ym=Cell.Centroid(2);
    xxo=Rfp.spotX; 
    yyo=Rfp.spotY; 
    xxt=Cfp.spotX; 
    yyt=Cfp.spotY;
    % 1) number of spots; enforce common ter pos
    N_ori=length(xxo);
    N_ter=length(xxt);  
    xt=mean(xxt);
    yt=mean(yyt);
    
    d_tm=((xt-xm).^2+(yt-ym).^2).^0.5; %distance ter to com
    d_ref=Cell.MinorAxisLength; %cell short axis as reference
    d_tm_Rel=d_tm/d_ref;
    
    % distances per ori or ter & typical cell 
    for ii=1:N_ori
        %ori to nearest ori (if N>1) & their shared 'pair center' 
        xo=xxo(ii);                  yo=yyo(ii); 
        if N_ori>1
            otheri=find([1:N_ori]~=ii);
            xxo_otheri=xxo(otheri); yyo_otheri=yyo(otheri);
            d_oo=((xxo_otheri-xo).^2+(yyo_otheri-yo).^2).^0.5; %distance to other ori
            [d_o2,nix]=min(d_oo);
            xp=xxo_otheri(nix); yp=yyo_otheri(nix);
            xo2=mean([xo xp]);  yo2=mean([yo yp]); %'pair center'
        else
            d_o2=0;
            xo2=xo; yo2=yo;
        end          
        %ori-ori-ter
        %pair center to nearest ter
        d_o2_t=((xo2-xt).^2+(yo2-yt).^2).^0.5; %distance pair center to ter
        %ori to com
        d_o2_m=((xo2-xm).^2+(yo2-ym).^2).^0.5; %distance pair center to com            
        %relative         
        d_o2_Rel=d_o2/d_ref;
        d_o2_t_Rel=d_o2_t/d_ref;
        d_o2_m_Rel=d_o2_m/d_ref;          
        ori_ori_ter_table(ii,:)=[ii d_tm_Rel d_o2_Rel d_o2_t_Rel d_o2_m_Rel] ;       
    end    