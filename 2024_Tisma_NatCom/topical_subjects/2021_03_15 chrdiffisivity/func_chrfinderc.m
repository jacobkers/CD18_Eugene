function chrout = func_chrfinderc(Ihup,chrnr2);
% output: label, area, center, left, right, length (all in microns), real chromosome number found
% This function uses background corrected hupA signal to find specified amount of chromosomes
% primarily used in chromosome line experiments with DnaCts
% Note20150925: here we also calculate the total intensity of the nucleoid
Ihup=double(Ihup);
Intbg=mean(mean(Ihup(:,1:2),2),1);
Ihup2=Ihup-Intbg;
if size(Ihup,2)>120;
    verti=max(Ihup,[],2);
    [~,mid]=max(verti,[],1);
    Ihupcrop=Ihup(max(mid-3,1):min(mid+3,25),:);
    hori=mean(Ihupcrop,1);
    [count,val]=hist(hori);
    [~,cmax]=max(count(1:floor(1/2*numel(count))));
    bg=val(cmax);
    Ihup=Ihup-bg;
    Hupthres=(max(hori)-bg)/4;
else
    Hupthres=max(Ihup(:))/4;% set threshold for binary image
end

Ibi=Ihup>Hupthres;
Ibi2=bwmorph(Ibi,'clean');
Ibi3=bwareaopen(Ibi2,64);
Ibi4=bwareaopen(Ibi2,100); % quality control, if not reached, the image is probably blank
if sum(Ibi4)==0;
    chrout=zeros(chrnr2,6);
    chrfound=0;
else
    [Lab,num]=bwlabel(Ibi3,4);
    % now pick out the size/center/left/right of each culster
    
    labparas=[];
    for i=1:num;
        [r,c]=find(Lab==i);
        area=numel(r);
        mid=mean(c); left=min(c); right=max(c);
        Totint=0;
        for xx=1:area;
            Totint=Totint+Ihup2(r(xx),c(xx));
        end
        labpara=[i area mid left right Totint Totint/area];
        labparas=cat(1, labparas, labpara);
    end
    labparas=-sortrows(-labparas,2); % line up area size at descending order
        if num==chrnr2;
            chrout=sortrows(labparas,3);% line up from left to right
            chrfound=chrnr2;
        elseif num<chrnr2;
            % add chromosomes from large to small for non-segregated chromosomes
            chrfound=num;
            chrout=labparas;
            add=chrnr2-num;
            multi=floor(add/num);
            if multi>0; 
                for j=1:multi; chrout=cat(1,chrout,labparas);
                end
            end
            addl=add-multi*num;
            if addl>0;
                chrout=cat(1,chrout,labparas(1:addl,:));
            end
            chrout=sortrows(chrout,3);
        else % now if too many clusters were found
            % first line up the biggest clusters as cores
            cores=labparas(1:chrnr2,:);
            chrfound=chrnr2;
            for k=chrnr2+1:num;
                leftdists=abs(labparas(k,4)-cores(:,5)); % the distance from left side of the cluster to the right side of all core clusters
                rightdists=abs(labparas(k,5)-cores(:,4)); % the opposite of the above;
                cldists=min([leftdists rightdists],[],2);
                [~,minloc]=min(cldists,[],1);
                newleft=min(labparas(k,4),cores(minloc,4));
                newright=max(labparas(k,5),cores(minloc,5));
                newmid=(newleft+newright)/2;
                newarea=labparas(k,2)+cores(minloc,2);
                cores(minloc,2:5)=[newarea newmid newleft newright];
            end
            chrout=sortrows(cores,3);
        end  
end
    chrout=cat(2,(1:1:chrnr2)',chrout(:,2)./12.5/12.5,chrout(:,3:5)./12.5,(chrout(:,5)-chrout(:,4))./12.5,ones(chrnr2,1).*chrfound,chrout(:,6:7));
end

