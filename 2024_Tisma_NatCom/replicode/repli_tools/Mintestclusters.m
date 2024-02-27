%test it
function Mintestclusters

 %% example analysis settings
    Psf=2.7;            %point spread function to use for peeling of spots
    Psf_meas=Psf;
    initval.Psf_est=Psf;
    ChipFract=1;        %fraction of peak height to subtract
    RelChangeStop=0.03; %stop criterio based on relative change of covered fraction
    sho=1;              %for plotting
    %% example blob picture
    psf0=2.5;
    hfz=25; blobno=20;
    blobposX=hfz/5*randn(blobno,1)+hfz;
    blobposY=hfz/5*randn(blobno,1)+hfz;
    chro_pic=zeros(2*hfz+1,2*hfz+1);
    for ii=1:blobno
        chro_pic=chro_pic+TwoDGaussNormPeak(chro_pic,blobposX(ii),blobposY(ii),psf0);
    end

AllSpotProps=PeelblobsFromImage(chro_pic,Psf_meas,0.98,0.03, 0);     
AllSpotProps=F006_CleanSpots(AllSpotProps,chro_pic,Psf_meas); 
ClusterBasics=Find_Clusters(AllSpotProps, initval);
[Clusters,ClusterPropsRow,ThisCellClusterTable,AllContours,ReconstructIm]=GetClusterDetails(ClusterBasics,chro_pic);

dum=1;

function Clusters=Find_Clusters(SpotProps, initval);
% 'Use this section for a Quicksheet'
%------------------------------------------------------------
    % sort spots by brightness; pick brightest one
    % find near ones for brightmost one from leftover list
    % for the ones found, keep finding new near ones until nochange; 
   
        % repeat for leftovers until no leftovers
        %JacobKers 2017 ----------------------------------------
        %AllSpotProps=[
        % 1 spotcount 
        % 2 Peak 
        % 3 Xpos 
        % 4 Ypos 
        % 5 Psf 
        %6 ThisSpotFraction 
        %7CoveredFraction 
        %8 RelChange]];
%------------------------------------------------------------[JK16]
% 'End of Quicksheet section'

 left_ones=SpotProps;
 Clusters=struct([]);
 N_clust=0;
 ClusterTable=[];
 while ~isempty(left_ones);    
    %build a new cluster
    N_clust=N_clust+1;
    Bri=left_ones(:,6); [Brisort,idx]=sort(Bri,'descend');
    thiscluster=left_ones(1,:);
    left_ones=left_ones(2:end,:);  %pick from stock
    [LC,~]=size(thiscluster);
    cluster_growth=1;
    while cluster_growth
        new_ones=[];
        for ii=1:LC  
            %for all existing elements of this cluster, 
            %,find near ones in the leftover stock
            x0=thiscluster(ii,3);
            y0=thiscluster(ii,4);
            leftx=left_ones(:,3);
            lefty=left_ones(:,4);
            rr=((leftx-x0).^2+(lefty-y0).^2).^0.5;
            SeparateWidth=2.5;
            near_ones_idx=find(rr<SeparateWidth*initval.Psf_est); %overlapping
            left_ones_idx=find(rr>=SeparateWidth*initval.Psf_est); 
          
            if ~isempty(left_ones_idx)
                new_ones=[new_ones; left_ones(near_ones_idx,:)]; %add 
                left_ones=left_ones(left_ones_idx,:);  %shrink remaining                 
            end
        end
        if ~isempty(new_ones) 
            thiscluster=[thiscluster; new_ones]; %add new ones to cluster 
        else
            cluster_growth=0;  %...or stop this cluster
        end 
    end
    %add values to cluster struct; go to next
    Clusters(N_clust).spotprops=thiscluster;    
 end
