function spot=Screen_spots(spot,pic,initval,Cell);
    %cleaning and screening spots (label-specific) image by image
    %remove (almost) merged ones
    %discard too blobby patterns or lack of black background
    %establish a 'cell area' testing pr spot::
        %max number nearby
        %if too high, remove if weak one
    if nargin<3
        initval.analyze_ori_ter=0;
        initval.analyze_replisome=1;
        pic=zeros(100,100);
        noiz=rand(100,100)*2;
        spot.spotX=[10:20:90  15:20:95 ceil(5*rand(1,6)+50)];
        spot.spotY=[15:20:95  15:20:95 ceil(5*rand(1,6)+50)];
        Lx=length(spot.spotX);
        spot.spot_ori_index=1:Lx;
        spot.spotContent=spot.spotX*0+2+8*rand(1,Lx);      
        for ii=1:Lx        
            pic(spot.spotY(ii),spot.spotX(ii))=spot.spotContent(ii);
        end
        pic=pic+noiz;
    end
    Lx=length(spot.spotX);
    spot.spot_ori_index=1:Lx;
    
    crit.mindistance=3; %if smaller, spots are considered one
    crit.pairdistance=15; %if mutually closest and smaller, spots are considered a pair
    crit.minpercentage=30/Cell.N_est; %minimum percentage spot of spot total
    %(adjusted for approximate cell number)
    crit.brightstars=30; %how much percentage intensity of whole should be 
                 %in good spots (to reject blobby images as a whole)
    %n_mx number per area
        
        
    %build screening recipe
    if initval.analyze_ori_ter        
    end
    if initval.analyze_replisome | initval.analyze_ori_ter       
    end   
    
    thrashit=0*spot.spotX;
    [xs,ys,Is,ix,orix,spot,Ns]=Update_spots(spot,thrashit); 
    
    

%% 1)first, remove identical or nearly-identical entries  
    ii=0;  stopit=0;
    
    if length(xs)~=0
    while~stopit       
        ii=ii+1;
        d_other=((xs-xs(ii)).^2+(ys-ys(ii)).^2).^0.5;
        nonself=(ix~=ii);
        Iii=Is(ii);
        doublecount=((d_other==0)&nonself);  
        nearlymerge=((d_other<=crit.mindistance)&nonself);
        otherweaker=((Is<=Iii)&nonself);
        to_be_removed=doublecount|(nearlymerge&otherweaker);
        if any(to_be_removed) %remove these 
            [xs,ys,Is,ix,orix,spot,Ns]=Update_spots(spot,to_be_removed); 
            ii=0;  %reset the counter
        end
        if ii==Ns, stopit=1; end
    end
    
    %reset the indices (thus, al former spots will be erased)
    spot.spot_ori_index=(1:length(spot.spot_ori_index));
    
 %% Now, don't shrink spots but screen for criteria
 %1 label strength of spots       
        Is_perc=Is/sum(Is)*100;
        for ii=1:Ns  %for each spot: 
            strong_enough=Is_perc(ii)>crit.minpercentage;
            spot.OK_strong(ii)=strong_enough;
        end
  
     
 %2) label pairs of strong spots
 spot=Find_strong_pairs(spot,crit);
 
    
  %3 evaluate 'blacknesss' of remaining picture
  %first, find all spots (including weak ones)
 
  switch 2
      case 1,goodspots=find(spot.OK_strong==1);   %take strong ones
      case 2,goodspots=find(spot.OK_strong*0+1);  %take all
  end
    goodcontent_perc=sum(spot.spotContent)/sum(pic(:))*100;
    if goodcontent_perc>crit.brightstars,
         spot.OK_nicepic=0*spot.OK_strong+1;; %accept picture
    else
        spot.OK_nicepic=0*spot.OK_strong; %reject picture
    end
    else
        spot=Nan_Fill(spot);
        
    end
    spot=orderfields(spot);  
    
     
    %plot menu; per type
    if nargin<3
        close all; pcolor(pic); shading flat, colormap hot; hold on
        %weakspots
        wk=find(spot.OK_strong==0);
        plot(spot.spotX(wk),spot.spotY(wk), 'yx', 'MarkerSize', 6); 
        %strongspots
        str=find(spot.OK_strong==1);
        plot(spot.spotX(str),spot.spotY(str), 'wo', 'MarkerSize', 8); 
        %build pairs
        sel=find(spot.OK_pair==1);
        LP=length(sel); pairX=zeros(LP,2); pairY=zeros(LP,2);
        for pp=1:LP
            id1=sel(pp); 
            id2=spot.OK_pair_idx(id1);
            pairX(pp,:)=[spot.spotX(id1) spot.spotX(id2)];
            pairY(pp,:)=[spot.spotY(id1) spot.spotY(id2)];           
        end       
        plot(pairX',pairY','wo-', 'MarkerSize',4,'MarkerFaceColor','w');
    end
    

    
 function spot=Find_strong_pairs(spot,crit);
 spot.OK_pair_idx=0*spot.spotX;
 spot.OK_pair=0*spot.spotX;
 spot.OK_pair_dist=0*spot.spotX; %mutually paired
 too_weak=~spot.OK_strong; %pick the strong ones:
 [str_x,str_y,str_I,str_ix,str_orix,~,Ns]=Update_spots(spot,too_weak); 
 str_N=length(str_x);
 if str_N>1;
        for ii=1:str_N  %for each spot: 
            %1) check possible pairs
            spot.OK_pair_idx(ii)=0;
            spot.OK_pair(ii)=0;
            [xp,yp,mindist, ip]=Find_nearest_neigbour(ii,str_x,str_y); 
            if mindist<=crit.pairdistance %possible pair
                [x2,y2,mindist, i2]=Find_nearest_neigbour(ip,str_x,str_y); 
                if i2==ii
                    idx1=str_orix(ii);  %original index spot 1
                    idx2=str_orix(ip);  %original index spot 2
                    spot.OK_pair(idx1)=1; %mutually paired
                    spot.OK_pair_idx(idx1)=idx2; %mutual pair
                    spot.OK_pair_dist(idx1)=mindist;
                end
            end
        end
 end   
dum=1; 
    
    
function [xs,ys,Is,ix,orix,outspot,Ns]=Update_spots(spot,thrashit);
            %return a selection of spots plus arrays for further working
            %note that repeated use will pass on original indices in
            %'spots' but not in array
            sel=find(thrashit==0);
            outspot.spotX=spot.spotX(sel);
            outspot.spotY=spot.spotY(sel);
            outspot.spotContent=spot.spotContent(sel);
            outspot.spot_ori_index=spot.spot_ori_index(sel);
            orix=spot.spot_ori_index(sel);
            
            xs=outspot.spotX;
            ys=outspot.spotY;
            Is=outspot.spotContent;
            Ns=length(xs); 
            ix=(1:Ns); 
            
  function spot=Nan_Fill(spot);
        spot.spotX=NaN; 
        spot.spotY=NaN;
        spot.spotContent=NaN;
        spot.spot_ori_index=NaN;
        spot.OK_strong=NaN;
        spot.OK_pair_idx=NaN;
        spot.OK_pair=NaN;
        spot.OK_pair_dist=NaN;
        spot.OK_nicepic=NaN;