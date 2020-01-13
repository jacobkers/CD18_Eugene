function StetsonCurve=prf_build_stetson_mask(Profile,edgepixels);
%JWJK_C*:-------------------------------------------------------------------
%Title: Build1DStetsonMask
%Project: various CD lab, Written by: Jacob
%Summary: %Build a hat-shaped mask, flat center and soft sloping edges, to
%suppress edges width a tunable width.
%Input: 1D-profile
%Output: masked 1D-profile 
%:JWJK_C*-------------------------------------------------------------------

%Build mask, see demo. JacobKers2018
if nargin<1
    Profile=(cos((1:30)/30*2*pi))';
end
%Build a flat-topped edge mask

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