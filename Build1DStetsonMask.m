function StetsonCurve=Build1DStetsonMask(Profile);
%Build mask, see demo. JacobKers2018
if nargin<1
    Profile=(cos((1:30)/30*2*pi))';
end
%Build a flat-topped edge mask
edgepixels=4;

LP=length(Profile);


%fill it with amplitude dependent masking
dropofflengthIn=edgepixels;
dropofflengthOut=edgepixels;

flatlength=LP-dropofflengthIn-dropofflengthOut;
EdgesIn=(hann(2*dropofflengthIn))';
EdgesOut=(hann(2*dropofflengthOut))';
StetsonCurve=[...
    EdgesIn(1:dropofflengthIn) ...
    ones(1,flatlength)...
    EdgesOut(dropofflengthOut+1:end)]'; 

if nargin<1
    close all;
    plot(Profile); hold on;
    plot(StetsonCurve,'b-'); hold on;
    plot(StetsonCurve.*Profile,'r-'); hold on;
    
end