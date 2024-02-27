function A001_Shell

user='Jacob';
batchindices=[1 2 3 4 5 6 7];
batchindices=[1];

for ii=1:length(batchindices)
    expno=batchindices(ii);
    switch user
    case 'Jacob'
        switch expno
            case 1, expi='BSG4595_1s_001';
            case 2, expi='BSG4595_1s_002';
            case 3, expi='BSG4595_1s_004';
            case 4, expi='BSG4595_10s_001_series1'; 
            case 5, expi='BSG4595_10s_001_series2';   
            case 6, expi='BSG4595_30s_001_series1'; 
            case 7, expi='BSG4595_30s_001_series2';   
        end
    end
initval=A000_ConfigExp(expi);

if 0, A010_GetSpots(initval); end
if 1, A020_BuildTraces(initval); end
end