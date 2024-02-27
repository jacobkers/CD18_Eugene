function loci=func_locifinder(Jspot1,sptnr,chrdata)
% the input image is already background subtracted
chrnr=size(chrdata,1);
            
%             Jspotbg=double(gaussf(Jspot1,10,'best'));
%             Jspot1=Jspot1-Jspotbg;
            Jspotbi=Jspot1>max(Jspot1(:))/3; % this is just a sensible thresholding, might need adjustment
            Jspotbi=bwmorph(Jspotbi,'clean');
            Jspotbi=bwareaopen(Jspotbi,8); % cluster at least with 8 pixels
            [Spotlabels,Spotnum]=bwlabel(Jspotbi,8);
            Spots=[];
            if Spotnum==0;
                loci=zeros(sptnr,4);
            else
                for m=1:Spotnum;              
                    [Spoty,Spotx]=find(Spotlabels==m);
                    Spotsum=0;Spotsumx=0;Spotsumy=0;
                    pixvals=[];
                    for m1=1:numel(Spoty);
                        pixval=Jspot1(Spoty(m1),Spotx(m1));
                        pixvals=cat(1,pixvals,pixval);
                        Spotsum=Spotsum+pixval;
                        Spotsumx=Spotsumx+pixval*Spotx(m1);
                        Spotsumy=Spotsumy+pixval*Spoty(m1);
                    end
                    Spotxc=Spotsumx/Spotsum;
                    Spotyc=Spotsumy/Spotsum;
                    pixsort=sort(pixvals);
                    pixmean=mean(pixsort(1:8));
                    Spots=cat(1,Spots,[Spotxc Spotyc pixmean m]);% coordinates and total intensity, in case too many were found                    
                end
                % if not enough spots found, then assume spots are not
                % segregated yet, so repeat them in the list.
                while size(Spots,1)<sptnr;
                    Spots=cat(1,Spots,Spots(1:min(size(Spots,1),sptnr-size(Spots,1)),:));
                end
                Spots=-sortrows(-Spots,3); % line up the brightness
                Spots=Spots(1:sptnr,:);
                Spots(:,4)=(1:1:sptnr)';
                % if no chromosome data, take that chromosomes cover the
                % whole cell
                if sum(chrdata(:))==0;
                    chrdata(:,4)=1;
                    chrdata(:,5)=size(Jspot1,2);
                end
                % now we match the Spots with the chromosomes and give a score for colocalization if loci are not inside but near one of them  
                loci=[];
                Spots(:,1:2)=Spots(:,1:2)./12.5;
                for j=1:sptnr;
                    scores=[];
                    for k=1:chrnr;
                        ld=Spots(j,1)-chrdata(k,4);rd=Spots(j,1)-chrdata(k,5);
                        if ld>=0 && rd<=0;
                            scores=cat(1,scores,1);
                        else scores=cat(1,scores,1-min(abs(ld),abs(rd))/(chrdata(k,5)-chrdata(k,4)));
                        end
                    end
                % now first assign the best matching loci
                    [chrS,chrN]=max(scores,[],1);
                    loci=cat(1,loci,[chrN Spots(j,1:2) chrS]);
                    loci=sortrows(loci,1);
                end
end