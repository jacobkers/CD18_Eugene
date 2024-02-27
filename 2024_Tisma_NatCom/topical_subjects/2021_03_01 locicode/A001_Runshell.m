function A001_Runshell
%this function lists the loci code programs in order of use

batch_set=[2.1];  %Timo '2019.12.5 - 2179 room temp local test
batch_set=[2.1 2.2 2.2 2.3 2.4 2.5];  %Timo '2019.12.5 - 2179 room temp remote test
%batch_set=2.2;  %Timo '2019.12.5 - 2179 room temp local test
batch_set=[2.1];  %Timo '2019.12.5 - 2179 room temp remote test
%batch_set=2.2;  %Timo '2019.12.5 - 2179 room temp local test
%% first, user-work
if 0
    for bi=1:length(batch_set)
        expstring=batch_set(bi);
        initval=A000_Lociicode_Get_PathsandExperiments(expstring);
        A005_ClickDriftMarkers(initval);
    end
end

%% then, auto-run
for bi=1:length(batch_set)
    expstring=batch_set(bi); 
    initval=A000_Lociicode_Get_PathsandExperiments(expstring);
    if 0, A010_FindDriftVector(initval);end
    if 0, A020_CorrectDrift(initval);end
    if 0, A040_Cellidentifier_LociData(initval);end
    if 0, A050_Cellcropper_Loci(initval);end
    if 0, A060_LociiTracking(initval);end
    if 0, A065_LociiTraceDynamics(initval);end
    if 1, A070_LociiShow(initval);end
end