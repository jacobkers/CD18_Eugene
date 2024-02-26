function A085_WF_BuildClusterDomainMaps(batchrunindex)
%JWJK_A:-------------------------------------------------------------------
%Two-dimensional (image-based) cluster analysis
%
%Summary: This function performs 2D cluster analysis on chromatin patterns of 
%cells and stores the result per cell, per cluster.
%
%Approach: we deconstruct the image by means of single-spot, psf-sized
%Gaussians. This is done by subtracting such Gaussians from the image by
%fist picking the brightest image point and subtracting one Gaussian of
%equal peak intensity there, then repeating this action on the residual
%image and continuing to do so until all intensity is covered by the subtracted 
%Gaussians.
%Next, the Gaussians are grouped in clusters. Each cluster consists of a
%small group of Gaussians(typically 1-5) or 'components' that are within one psf 
%distance of one of the others, such that the cluster forms an optically connected shape. 
%Likewise, all components of one cluster are at least one psf away form
%those of another.
%These clusters are then parametrized by their center-of mass position,
%content, radius of gyration and so on. Finally, their angular position with respect to  
% the chromatin geometrical center (as was determined via the 1D density curve analysis)
% is determined and stored.
%
%Input: data in .mat files stored in former analysis steps.
%
%Output: Data is presented as scatter plots and histograms. Saved are 
%Tabular excel files, .mat and summary plots.
%
%
%:JWJK_A-------------------------------------------------------------------


close all;
if nargin<1,batchrunindex=11.1;end
sho=0;
initval=A000__WF_Get_JacobPathsandExperiments(batchrunindex);

imoutdir=strcat(initval.resultpath,'A085_DomainMaps',initval.DirSep);
if isdir(imoutdir), rmdir(imoutdir,'s');  end
mkdir(imoutdir);
mkdir(strcat(imoutdir,initval.DirSep,'Tiffs'));
if sho
    mkdir(strcat(imoutdir,initval.DirSep,'Plots'));
end


allframes=length(initval.Cell_Labels);


goodcount=0;
MatFilePath=strcat(initval.resultpath,'ResultsPerCellMatlab',initval.DirSep);
for jj=1:allframes    
    cellno=char(initval.Cell_Labels{jj});     
    CellName=strcat('ResultsOfCell',cellno,'.mat'); 
        load(strcat(MatFilePath,CellName));   
        if GeneralCellProps.Okayproduct|initval.PassAllCells
        %disp('good cell')
        disp(strcat(CellName,'ClusterAnalysis..', num2str(allframes-jj+1), 'cells to go'));
        goodcount=goodcount+1;        
        cellno=char(initval.Cell_Labels{jj});    
        NumCellLabel=BuildNumericCellLabel(cellno);
        CellSpecs=[goodcount NumCellLabel];      
        [~,CL]=size(Clusters);
        
        
        %if strcmp(cellno,'100285'),sho=1;, else sho=0;end
            close all;
            if sho
            subplot(2,2,1);
            pcolor(chro_pic); shading flat; colormap bone; hold on;
            title(CellName);
            end
            clusterDist=zeros(CL,1);
            clusterBP=zeros(CL,1);
            for cc=1:CL %for all clusters               
                %spotProps=
                    % [spotcount 
                    % Peak 
                    % Xpos 
                    % Ypos 
                    % Psf 
                    % ThisSpotFraction(spotcount) 
                    % CoveredFraction(spotcount) 
                    % RelChange]]; 
                if sho
                    subplot(2,2,1);
                    plot(Clusters(cc).COM_X,Clusters(cc).COM_Y,'ro'); hold on;
                    plot(Clusters(cc).spotprops(:,3),Clusters(cc).spotprops(:,4),'wo-');
                end
                clusterDist(cc)=Clusters(cc).COM_DistPerc;
                clusterBP(cc)=Clusters(cc).C_perc;       
            end
            
            [clusterDist,idx]=sort(clusterDist);
            clusterBP=clusterBP(idx);
            clusterCumBP=100*cumsum(clusterBP);
            clusterCumBP=clusterCumBP/max(clusterCumBP)*100;  %small correction
                        
            BPipolax=1:100;
            DistanceRe_IPol=interp1(clusterCumBP,clusterDist,BPipolax);
            %build the heat map per image
            DistanceHeatMap=zeros(100,100); 
            clustcount=1; dist=0;
            for ii=1:100
                for jj=1:100
                    [val,nearest_i_clust_idx]=min(abs(ii-clusterCumBP));
                    dist_i=clusterDist(nearest_i_clust_idx);
                    
                    [val,nearest_j_clust_idx]=min(abs(jj-clusterCumBP));
                    dist_j=clusterDist(nearest_j_clust_idx);
                    deltadist=abs(dist_i-dist_j);
                    DistanceHeatMap(100-jj+1,ii)=abs(100-deltadist);
                end
            end 
            if sho
                subplot(2,2,2);
                pcolor(DistanceHeatMap); shading flat; colormap bone; hold on;
                pause(0.01);
                %[~]=ginput(1);  
            end
            OutName=strcat('DomainMap',CellName(1:end-4));
        
            %build a tiff picture
            imout=uint16(10*(flipud(DistanceHeatMap)-1));
            imwrite(imout,strcat(imoutdir,'Tiffs',initval.DirSep,OutName,'.tif'),'tif');
            if sho
                saveas(gcf,strcat(imoutdir,'Plots',initval.DirSep, OutName,'.jpg')); 
                pause(0.01);
            end    
  
    end
end

              
  
  
        
        