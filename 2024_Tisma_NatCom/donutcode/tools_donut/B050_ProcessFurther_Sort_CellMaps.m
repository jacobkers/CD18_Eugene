function [ChromosomeProfilesNw,OutSortIdx]=G00_ProcessFurther_Sort_CellMaps(ChromosomeProfiles,sortstyle,impose_idx);
    %This function counts the minima per row and sorts the array in that
    %order
    %1) clean
    sel=find(isnan(ChromosomeProfiles));
    ChromosomeProfiles(sel)=0;
    
    [rr,cc]=size(ChromosomeProfiles);
    
    if nargin<3, 
        impose_idx=1:rr; 
    end
    
    switch sortstyle
        case 'None'  %keep as-is
           ChromosomeProfilesNw=ChromosomeProfiles;
           OutSortIdx=1:rr;
        case 'Imposed'  %use former sorting index
            ChromosomeProfilesNw=ChromosomeProfiles(impose_idx,:);  
            OutSortIdx=impose_idx;
        case 'Std'  %sort by variation in profile
           stdcount=zeros(1,rr);
            for ii=1:rr  
                prf=ChromosomeProfiles(ii,:);
                stdcount(ii)=nanstd(prf);            
            end
            [~,stdSortIdx]=sort(stdcount,'ascend');
            ChromosomeProfilesNw=ChromosomeProfiles(stdSortIdx,:);   
            OutSortIdx=stdSortIdx;
        case 'Peak' %sort by peak value per profile
           pkcount=zeros(1,rr);
            for ii=1:rr  
                prf=ChromosomeProfiles(ii,:);
                pkcount(ii)=nanmax(prf);            
            end
            [~,pkSortIdx]=sort(pkcount,'ascend');
            ChromosomeProfilesNw=ChromosomeProfiles(pkSortIdx,:);   
            OutSortIdx=pkSortIdx;
         case 'PeakPos' %sort by peak value per profile
           pk_pos=zeros(1,rr);
            for ii=1:rr  
                prf=ChromosomeProfiles(ii,:);
                [~,idx]=nanmax(prf); 
                pk_pos(ii)=idx          
            end
            [~,pkSortIdx]=sort(pk_pos,'descend');
            ChromosomeProfilesNw=ChromosomeProfiles(pkSortIdx,:);   
            OutSortIdx=pkSortIdx;
        case 'Minima'  %sort by number of minima
            %1) count minima
           mincount=zeros(1,rr);
            for ii=1:rr  
                prf=ChromosomeProfiles(ii,:);
                prf_left=prf(1:end-2);
                prf_mid=prf(2:end-1);
                prf_right=prf(3:end);
                mincount(ii)=nansum(1.0*((prf_mid<prf_left)&(prf_mid<prf_right)));            
            end
            [~,MinimaSortIdx]=sort(mincount,'ascend');
             ChromosomeProfilesNw=ChromosomeProfiles(MinimaSortIdx,:);
            %2) normalize     
             for ii=1:rr
                 prf=ChromosomeProfilesNw(ii,:);
                 prf_n=prf/nanstd(prf)*100;
                 ChromosomeProfilesNw(ii,:)=prf_n;
             end
           OutSortIdx=MinimaSortIdx;
        case 'Similarity'
           %this will align based on the xxxth quarter (i.e., around the
           %start label)
           [ChromosomeProfilesNw, OutSortIdx]=group_rows_by_similarity(ChromosomeProfiles);           
    end
    
    dum=1;