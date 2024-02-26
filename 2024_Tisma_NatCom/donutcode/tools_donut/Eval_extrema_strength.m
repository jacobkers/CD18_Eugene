function  DistancePerc=Eval_extrema_strength(extrema,DensityAxis,DensityCurve,perc) ;
    %wat range around this extremum gives us 'perc' content?
    tripledDensityCurve=repmat(DensityCurve,1,3);
    contentperpoint=nanmax(DensityAxis)/length(DensityAxis);
    shft=length(DensityCurve);
    idxes=extrema.Idx+shft;
    LE=length(idxes);
    DistancePerc=zeros(1,0);
    for ii=1:LE
        halfdist=0; content=0;
        idx=idxes(ii);
        if ~isnan(idx)
            while content<perc
                halfdist=halfdist+1;
                content=nansum(tripledDensityCurve(idx-halfdist:idx+halfdist));
            end
           DistancePerc(ii)=(2*halfdist)*(contentperpoint);
           if DistancePerc(ii)<=1  %single-point approach
           DistancePerc(ii)=tripledDensityCurve(idx)*contentperpoint/perc;
%            else
%                DistancePerc(ii)=NaN;
           end
       end
    end
    dum=1;
