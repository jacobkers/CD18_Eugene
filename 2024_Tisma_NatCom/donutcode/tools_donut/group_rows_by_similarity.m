function [outmap, OutSortIdx]=group_rows_by_similarity(inmap)
           outmap=inmap;
           [rr,cc]=size(inmap);
           SortIdx=1:rr;

           %first, pick first one:
           preselect=zeros(1,rr);
            for ii=1:rr  
                prf=inmap(ii,:);
                preselect(ii)=nanstd(prf);            
            end
            [preselect,preSortIdx]=sort(preselect,'descend');
            ChromosomeProfilesPre=inmap(preSortIdx,:);   
           PickedProfiles=ChromosomeProfilesPre(1,:);
           
           %then, find similar ones  
           PickedOriIdx=[];
           RemainingProfiles=outmap(2:end,:);
           RemainingOriIdx=SortIdx;
           lookinsection=cc; %  lookinsection=ceil(cc/4);
           for ii=1:rr -1
               RefProfile=PickedProfiles(end,:);  %last picked one
               [rr2,~]=size(RemainingProfiles);
               DifCurve=sum(((RemainingProfiles(:,1:lookinsection)-repmat(RefProfile(1:lookinsection),rr2,1)).^2),2);
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
           outmap=PickedProfiles;
           OutSortIdx=PickedOriIdx;