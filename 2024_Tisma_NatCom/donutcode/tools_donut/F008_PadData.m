function [axz_pad,datamap_pad]=F008_PadData(axz,datamap,initval,modus)
     %This function add repeats to front and back end of data. 'map' can be
     %empty or 1D, modus: x-axis can be repeated, also Y-axis (assuming
     %monotonic increase)
     [rr,cc]=size(datamap);
     pd=ceil(initval.Padcurves/100*cc);
     maxax=max(axz);
     axz_pad=[axz(end-pd+1:end)-maxax axz(2:end-1) axz(1:pd)+maxax];
     if ~isempty(datamap)
        if strcmp(modus,'Xrepeat')
            datamap_pad=[datamap(:,end-pd+1:end) datamap(:,2:end-1) datamap(:,1:pd)]; end
        if strcmp(modus,'XYrepeat')
            Y_up=repmat(datamap(:,end),1,pd);
            datamap_pad=[datamap(:,end-pd+1:end)-Y_up datamap(:,2:end-1) datamap(:,1:pd)+Y_up]; end
     else
         datamap_pad=[];
     end