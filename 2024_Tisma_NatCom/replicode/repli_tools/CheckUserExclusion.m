
function notexcluded=CheckUserExclusion(CellName, initval);
notexcluded=1;
if isfield (initval, 'exclusionlist')
sel=find(strcmp(initval.exclusionlist,CellName));
if ~isempty(sel)
    notexcluded=0;
end
end