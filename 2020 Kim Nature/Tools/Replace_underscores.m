function outstring=Replace_underscores(instring)
if nargin<1
    instring='jacob_test_1'
end
outstring=instring;
idx=strfind(instring, '_');
if~isempty(idx)
    outstring(idx)='-';
end
