function A110_AnalyzeClusterGeometries(batchrunindex)

close all;
options.aligncellcenters=1;

if nargin<1,batchrunindex=4;end
initval=A000_Repli_Init(batchrunindex);
disp(initval.expi);
MatFilePath=strcat(initval.resultpath,'ResultsPerCellMatlab',initval.DirSep);

load(strcat(MatFilePath,strcat('__MovieList.mat')));
[~,LC]=size(MovieList);

for ii=1:LC 
    CellName=char(MovieList(ii).CellNames); disp(CellName);
    LF=length(MovieList(ii).CellFrames);     %per frame
    COM=[]; TER=[];
    for jj=1:LF        
        if jj>1, Cfp_former=Cfp; Rfp_former=Rfp; Cell_former=Cell; end  %see also line 30
        CellFrameName=char(MovieList(ii).CellFrames(jj)); disp(CellFrameName);
        %info files earlier analysis
        CellInfoPath=strcat(MatFilePath,'ResultsOfCell',CellFrameName);
        load(strcat(CellInfoPath,'_Spots.mat'), 'Cfp','Rfp','Cfp_pic','Rfp_pic');
        load(strcat(CellInfoPath,'_Clusters.mat'),'Clusters','chro_pic','celledge_pic');
        load(strcat(CellInfoPath,'_Cellshape.mat'),'Cell');        
        load(strcat(CellInfoPath,'_SpotPairs.mat'), 'oridist_table','oripair_table');        
        
            %oripair_table contains for this frame:
            %[ii [xo xp xo2] : 1234 xpos  of ori spot, neighbour, common center
            %    [yo yp yo2] : 567 ypos  of ori spot, neighbour, common center
            %     xm ym]     : 89 cell center-of-mass for reference
            %     xf yf      : 1011 xpos, ypos  of nearest ori spot former frame
        
        LP=size(oripair_table);
        [~,LC]=size(Clusters);
        for ci=1:LC   %get clusters
            clusters_x(ci)=Clusters(ci).COM_X-Cell.Centroid(1);
            clusters_y(ci)=Clusters(ci).COM_Y-Cell.Centroid(2);
            clusters_I(ci)=Clusters(ci).C_perc;
            
        end
        mrkr=ceil(30*clusters_I);
        for kk=1:LP   %For each ori pair
            xo=oripair_table(kk,2); yo=oripair_table(kk,5);
            xp=oripair_table(kk,3); yp=oripair_table(kk,6);
            %Get center&angle
            xo2=oripair_table(kk,4); yo2=oripair_table(kk,7);
            pair_angle=(atan2(yp-yo,xp-xo))/pi*180;
            
           
            %rotate ori pair
            [xxr,yyr]=Rotate_Points(xo2,yo2,[xo xp],[yo yp],-pair_angle);
            %rotate clusters
            [cxr,cyr]=Rotate_Points(xo2,yo2,clusters_x,clusters_y,-pair_angle);
            
            pair_sep=xxr(2)-xxr(1);
            %plot([xo xo2 xp],[yo yo2 yp],'ko-'); axis equal; hold on;
            %plot([xo],[yo],'ro-', 'MarkerFaceColor', 'r');
            if pair_angle>0 %use only half of the two possibilities
                plot(xxr-xo2, pair_sep, 'ro'); hold on;
                for ci=1:LC 
                    %plot(cxr(ci)-xo2, cyr(ci)-yo2, 'ko','MarkerSize',mrkr(ci)); hold on;
                    plot(cxr(ci)-xo2, pair_sep, 'ko','MarkerSize',mrkr(ci)); hold on;
                end
                xlabel('ori-bar projected position');
                ylabel('ori-pair separation');
            end
            %Get cluster coords
            %Rotate coords; calculate 'ori bar coordinate'
            %collect or scatter in plot
            %[~]=ginput(1);
        end
    end
end
