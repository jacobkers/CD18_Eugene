function mxi=JKD1_PRF_get1Dpeakssimple(data,sig)
%Really simple peakfinding routine
%Find peaks in a profile 'data'. A peak is a local maximum that is higher
%than a certain treshold.

%input: 1D data and a treshold expressed in sigmas of the standard
%deviation
%output: indices of the accepted peaks
%Jacob Kers 2013
%--------------------------------------------------------------------------

if nargin<2  %For testing puroposes
sig=0.5;
data=JK00_DEMODATA_Peaks;
end

md=median(data);
st=std(data);

all=find(data(2:end-1)>data(1:end-2) &data(2:end-1)>data(3:end))+1;  %indices of local maxes 1D
pks=data(all);
sel=find((pks>md+sig*st));
mxi=all(sel);

if nargin<2
plot(data,'o-'); hold on;
plot(mxi, data(mxi),'ro', 'MarkerFaceColor', 'r');
title('''F110 Get1DpeaksSimple'' DEMO')
end

