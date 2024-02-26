function firstpeaks=prf_get_peaks_by_fixedtreshold(prf,tres);
%super simple peak detection; we assume peaks are pre-selected, flat bottom
%and not too different

if nargin<2  %DEMO
    close all;
    prf=prf_make_demo_curves('multipeaks');
    prf=smooth(prf,3)-min(prf);
    tres=0.5*max(prf);
end



localmaxidx=find(prf(2:end-1)>prf(1:end-2) &prf(2:end-1)>prf(3:end))+1;  %local maxes 1D
sel=find(prf(localmaxidx)>tres);
firstpeaks=localmaxidx(sel);

if nargin<2 %optional plot menu-------------------------------------------
    close all;
    plot(prf, 'o-'); hold on;
    stem(firstpeaks, prf(firstpeaks), 'ko', 'MarkerFace', 'k'); %vertical lines at peak positions
    [~]=ginput(1); close(gcf);
end

