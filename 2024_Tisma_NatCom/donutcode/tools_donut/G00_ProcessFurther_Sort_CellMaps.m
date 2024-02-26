function [ChromosomeProfilesNw,OutSortIdx]=G00_ProcessFurther_Sort_CellMaps(ChromosomeProfiles,sortstyle);
    %This function counts the minima per row and sorts the array in that
    %order
    %1) clean
    sel=find(isnan(ChromosomeProfiles));
    ChromosomeProfiles(sel)=0;
    
    [rr,cc]=size(ChromosomeProfiles);
    
    switch sortstyle
        case 'None';
           ChromosomeProfilesNw=ChromosomeProfiles;
           OutSortIdx=1:rr;
        case 'Minima',
           mincount=zeros(1,rr);
     %first, sort minima per frame
            for ii=1:rr  
                prf=ChromosomeProfiles(ii,:);
                prf_left=prf(1:end-2);
                prf_mid=prf(2:end-1);
                prf_right=prf(3:end);
                mincount(ii)=nansum(1.0*((prf_mid<prf_left)&(prf_mid<prf_right)));            
            end
            [mincount,MinimaSortIdx]=sort(mincount,'ascend');
             ChromosomeProfilesNw=ChromosomeProfiles(MinimaSortIdx,:);
       %second, normalize     
             for ii=1:rr
                 prf=ChromosomeProfilesNw(ii,:);
                 prf_n=prf/nanstd(prf)*100;
                 ChromosomeProfilesNw(ii,:)=prf_n;
             end
           OutSortIdx=MinimaSortIdx;
        case 'Similarity' 
           figure;
           ChromosomeProfilesNw=ChromosomeProfiles;
           SortIdx=1:rr;
           PickedProfiles=ChromosomeProfilesNw(1,:);
           PickedOriIdx=[];
           RemainingProfiles=ChromosomeProfilesNw(2:end,:);
           RemainingOriIdx=SortIdx;
           for ii=1:rr -1
               RefProfile=PickedProfiles(end,:);  %last picked one
               [rr2,~]=size(RemainingProfiles);
               DifCurve=sum(((RemainingProfiles-repmat(RefProfile,rr2,1)).^2),2);
               [~, MostSimilarIdx]=min(DifCurve);

               MostSimilarProfile=RemainingProfiles(MostSimilarIdx,:);
               MostSimilarOriIndex=RemainingOriIdx(MostSimilarIdx);
               
               OtherIdx=find((1:rr2)~=MostSimilarIdx);
               OtherOriIndex=RemainingOriIdx(OtherIdx);
               PickedProfiles=[PickedProfiles; MostSimilarProfile]; 
               RemainingProfiles=RemainingProfiles(OtherIdx,:); 
               
               RemainingOriIdx=RemainingOriIdx(OtherIdx);
               PickedOriIdx=[PickedOriIdx MostSimilarOriIndex];        
           end
           ChromosomeProfilesNw=PickedProfiles;
           OutSortIdx=PickedOriIdx;
    end
    
    dum=1;