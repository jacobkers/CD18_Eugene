function stepsout=Merge_Edges(stepsin,kymo_stepped,kymo);
    %this function merges step-trains related to 'soft steps' per row by performing a
    %weighted merge of step coordinates (a cheap version of bluring the
    %locations and finding a peak....
    %per line: merge steps if closer than psf; 
    %get new locations
    %input: [stepidx,posx posy levelbefore levelafter];
    if nargin<1
        outname=strcat(pwd,'\mergetest.mat');
        load(outname,'kymo','kymo_stepped','stepsin');
        subplot(1,2,1); pcolor(kymo); shading flat;
        subplot(1,2,2); pcolor(kymo_stepped'); shading flat
        colormap hot;
    end
    [rr,cc]=size(kymo);
    for ii=1:rr
        %find asociated locations
        %per step; check near ones; 
        %merge
    end
    stepsout=stepsin;
