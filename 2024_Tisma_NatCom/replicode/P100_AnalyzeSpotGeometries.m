function P100_AnalyzeSpotGeometries(batchrunindex)
%This function analyzes relations beteen points and frames to follow spots
%and pairs of spots. We can do so with ori and ter (rfp/cfp) and with
%replisome (rfp only).
close all;
options.aligncellcenters=1;

initval=A000_Repli_Init(batchrunindex);
disp(initval.expi);
MatFilePath=strcat(initval.resultpath,'ResultsPerCellMatlab',initval.DirSep);
plotoutdir=strcat(initval.resultpath,'SpotAnalysis',initval.DirSep);  
if isdir(plotoutdir), 
    rmdir(plotoutdir,'s');  
end
mkdir(plotoutdir);

%% analysis section
load(strcat(MatFilePath,strcat('__MovieList.mat')));
[~,LC]=size(MovieList);
close all;
for ii=1:LC     
    CellName=char(MovieList(ii).CellNames);  
    ThisMovie=MovieList(ii);               
    if initval.analyze_ori_ter
        Analyze_ori_ter_relations(MatFilePath,ThisMovie,options);
    end
    if initval.analyze_replisome
        Analyze_replisome_relations(MatFilePath,ThisMovie,options);
    end
end


%% plot menu development
close all
for ii=1:LC 
    CellName=char(MovieList(ii).CellNames);  
    ThisMovie=MovieList(ii);     
    LF=length(ThisMovie.CellFrames);     %per frame
    for jj=1:LF 
        CellFrameName=char(MovieList(ii).CellFrames(jj));
        CellInfoPath=strcat(MatFilePath,'ResultsOfCell',CellFrameName);
        if initval.analyze_ori_ter                    
            Plot_ori_ter_general_panels(MatFilePath,CellFrameName);                  
        end 
        if initval.analyze_replisome                    
            Plot_replisome_panels(MatFilePath,CellFrameName,ii);                  
        end 
        
        
    end
        %image per cell
        if initval.analyze_ori_ter       
            figure(ii); subplot(2,2,1);   
            plot(COM(:,1),COM(:,2),'k-');  %cell center off mass
            plot(TER(:,1),TER(:,2),'b-');  %route of brightest spot    
            saveas(gcf,strcat(plotoutdir, CellFrameName,'_Ori_Ter_Dance.jpg')); 
        end
    %close(gcf);
    pause(0.1); %[~]=ginput(1);
end


function [ori_pairsX,ori_pairsY,FormerNearXYCoords]=Form_LookingBackPairs(Rfp,Rfp_former,Cell,Cell_former,options);
    %JWJK_B:-------------------------------------------------------------------
    %Title: Form time-pairs
    %Summary: %this function forms coordinate pairs of a spot and its nearest spot in the
    %former image; it is assumed these represent the same object in time     
    %''Rfp'' with fields: spotY, spotX,spotContent, same for fomrer image    
    %Output: separate X and Y coordinates of associated pairs and an xy
    %list for future use
    %:JWJK_B-------------------------------------------------------------------
    ori_pairsX=[];ori_pairsY=[];
    xm=Cell.Centroid(1);
    ym=Cell.Centroid(2);
    xmf=Cell_former.Centroid(1);
    ymf=Cell_former.Centroid(2);
    
    %drop the bad spots
    goodspots=find(Rfp.spotOk);
    Rfp.spotX=Rfp.spotX(goodspots);
    Rfp.spotY=Rfp.spotY(goodspots);
    Rfp.spotContent=Rfp.spotContent(goodspots);
    Rfp.spotOk=Rfp.spotOk(goodspots);
    
    %drop the bad spots
    goodspots=find(Rfp_former.spotOk);
    Rfp_former.spotX=Rfp_former.spotX(goodspots);
    Rfp_former.spotY=Rfp_former.spotY(goodspots);
    Rfp_former.spotContent=Rfp_former.spotContent(goodspots);
    Rfp_former.spotOk=Rfp_former.spotOk(goodspots);
    
    if ~options.aligncellcenters
    xf=Rfp_former.spotX;        yf=Rfp_former.spotY;
    x=Rfp.spotX;                y=Rfp.spotY;
    else
        xf=Rfp_former.spotX-xmf;        yf=Rfp_former.spotY-ymf;
        x=Rfp.spotX-xm;                 y=Rfp.spotY-ym;
    end
    Lp=length(y);
    FormerNearXYCoords=[];
    for ii=1:Lp                          %for every new spot
        xi=x(ii);yi=y(ii);
        rr=((xf-xi).^2+(yf-yi).^2).^0.5; %distance to all former
        [mindist,nearest]=min(rr);    
        ori_pairsX=[ori_pairsX [xi xf(nearest)]'];
        ori_pairsY=[ori_pairsY [yi yf(nearest)]'];
        FormerNearXYCoords(ii,:)=[xf(nearest) yf(nearest) mindist];
    end

 function  Analyze_ori_ter_relations(MatFilePath,ThisMovie,options);
        %oripair_table contains for this frame:
            %[ii [xo xp xo2] : 1234 xpos  of ori spot, neighbour, common center
            %    [yo yp yo2] : 567 ypos  of ori spot, neighbour, common center
            %     xm ym]     : 89 cell center-of-mass for reference
            %     d_min      :distance between label pair
            %     xf yf      : 1011 xpos, ypos  of nearest ori spot former frame
            %     d_min_former: distance to this former spot
             LF=length(ThisMovie.CellFrames);     %per frame
            for frameno=1:LF      
                CellFrameName=char(ThisMovie.CellFrames(frameno));
                CellInfoPath=strcat(MatFilePath,'ResultsOfCell',CellFrameName); 
                load(strcat(CellInfoPath,'_Clusters.mat'),'Clusters','chro_pic','celledge_pic');
                load(strcat(CellInfoPath,'_Cellshape.mat'),'Cell');    
                load(strcat(CellInfoPath,'_Spots.mat'), 'All_labels');
                Cfp=All_labels.Cfp;
                Cfp_pic=All_labels.Cfp_pic;
                Rfp=All_labels.Rfp;
                Rfp_pic=All_labels.Rfp_pic;
                if frameno==1, 
                    Cfp_former=Cfp; 
                    Rfp_former=Rfp; 
                    Cell_former=Cell; 
                end %just fillers        
                oridist_table=Get_ori_ori_ter_distances(Cfp,Rfp,Cell);  
                oripair_table=Get_labelpair_distances(Rfp,Cell,options);    
                [~,~,ori_formernearxycoords]=Form_LookingBackPairs(Rfp,Rfp_former,Cell,Cell_former,options); 
                oripair_table=[oripair_table ori_formernearxycoords];
                save(strcat(CellInfoPath,'_SpotPairs.mat'), 'oridist_table','oripair_table'); 
            end
            
  function  Analyze_replisome_relations(MatFilePath,ThisMovie,options);
                %replipair_table contains for this frame:
                %[ii [xo xp xo2] : 1234 xpos  of ori spot, neighbour, common center
                %    [yo yp yo2] : 567 ypos  of ori spot, neighbour, common center
                %     xm ym]     : 89 cell center-of-mass for reference
                %     xf yf      : 1011 xpos, ypos  of nearest ori spot former frame
                 LF=length(ThisMovie.CellFrames);     %per frame
                for frameno=1:LF      
                    CellFrameName=char(ThisMovie.CellFrames(frameno));
                    CellInfoPath=strcat(MatFilePath,'ResultsOfCell',CellFrameName); 
                    load(strcat(CellInfoPath,'_Clusters.mat'),'Clusters','chro_pic','celledge_pic');
                    load(strcat(CellInfoPath,'_Cellshape.mat'),'Cell');    
                    load(strcat(CellInfoPath,'_Spots.mat'), 'All_labels');
                    Rfp=All_labels.Rfp;
                    
                    
                    if frameno==1, 
                        Rfp_former=Rfp; 
                        Cell_former=Cell; 
                    end %just fillers        
                    replipair_table=Get_labelpair_distances(Rfp,Cell,options); 
                    [~,~,ori_formernearxycoords]=Form_LookingBackPairs(Rfp,Rfp_former,Cell,Cell_former,options); 
                    replipair_table=[replipair_table ori_formernearxycoords];
                    %screen pairs; add to full pair-table
                    %are the two spots mutually closest?
                    %is their distance small enough OR is their midpoint
                    %distance close enough to a former valid pair?
                    save(strcat(CellInfoPath,'_SpotPairs.mat'),'replipair_table'); 
                end
                %add here some info on movie: timelines of spots
                
                
                
                
                
     function Plot_ori_ter_general_panels(MatFilePath,CellFrameName,figureno);
            %plot ter-ori relations
            CellInfoPath=strcat(MatFilePath,'ResultsOfCell',CellFrameName);
            load(strcat(CellInfoPath,'_SpotPairs.mat'), 'oridist_table','oripair_table');
            load(strcat(CellInfoPath,'_Clusters.mat'),'Clusters','chro_pic','celledge_pic');
            load(strcat(CellInfoPath,'_Cellshape.mat'),'Cell');    
            load(strcat(CellInfoPath,'_Spots.mat'), 'All_labels');
            Cfp=All_labels.Cfp;
                Cfp_pic=All_labels.Cfp_pic;
                Rfp=All_labels.Rfp;
                Rfp_pic=All_labels.Rfp_pic;
            %% test representation actions        
               figure(figureno); subplot(2,2,1);
               timepair_x=[oripair_table(:,2) oripair_table(:,10)]';
               timepair_y=[oripair_table(:,5) oripair_table(:,11)]';
               cellcomx=Cell.Centroid(1);
               cellcomy=Cell.Centroid(2);
               
               if ~options.aligncellcenters
                    COM(jj,:)=[cellcomx cellcomy];
                    TER(jj,:)=[Cfp.spotX(1),Cfp.spotY(1)];
               else
                    COM(jj,:)=[0 0];
                    TER(jj,:)=[Cfp.spotX(1)-Cell.Centroid(1),Cfp.spotY(1)-Cell.Centroid(2)];
               end
               plot(timepair_x,timepair_y, 'r-'); hold on;            
                if jj==LF  %plot last contour
                    [rr,cc]=size(celledge_pic); [XX,YY]=meshgrid(1:cc,1:rr);
                    sel=find(celledge_pic==1);
                    Cx=XX(sel); Cy=YY(sel);  
                    if options.aligncellcenters
                        Cx=Cx-cellcomx;
                        Cy=Cy-cellcomy;
                    end
                end

                if jj==LF  %plot last contour 
                    figure(figureno); subplot(2,2,1);
                    plot(Cx,Cy,'ko','MarkerSize',3, 'MarkerFaceColor','y'); hold on; 
                    axis equal;
                    figure(figureno); subplot(2,2,2);
                    plot(Cx,Cy,'ko','MarkerSize',3, 'MarkerFaceColor','y');  hold on;
                    axis equal;
                end

                figure(figureno); subplot(2,2,1);
                plot(TER(jj,1),TER(jj,2),'bo-','MarkerSize',sqrt(jj+1)+1,'MarkerFaceColor','b');  %route of brightest spot
                plot(COM(jj,1),COM(jj,2),'ko-','MarkerSize',sqrt(jj+1)+1,'MarkerFaceColor','k');  %route of brightest spot
                plot(oripair_table(:,2),oripair_table(:,5),'ro','MarkerSize',sqrt(jj+1)+1,'MarkerFaceColor','r');
                title(strcat(initval.expi,'-',CellName)); 
                % legend('cellcom'   ,      'ter'     ,    'ori'  ,...
                %    'Location','SouthOutside');      
                figure(figureno); 
                subplot(2,2,2);        
                    lp=length(oripair_table(:,1));
                    oo_pairx=[oripair_table(:,2:3)]'; xo=oripair_table(:,4);
                    oo_pairy=[oripair_table(:,5:6)]'; yo=oripair_table(:,7);
                    om_pairx=[oripair_table(:,3) oripair_table(:,8)]'; xm=oripair_table(:,8);
                    om_pairy=[oripair_table(:,6) oripair_table(:,9)]'; ym=oripair_table(:,9);
                    plot(oo_pairx,oo_pairy,'ro-','MarkerSize',sqrt(jj+1)+1, 'MarkerFaceColor','r'); hold on; 
                    plot(xo,yo,'ro','MarkerSize',sqrt(jj+1)+1, 'MarkerFaceColor','w'); hold on; 
                    %plot(om_pairx,om_pairy,'k-','MarkerSize',sqrt(jj+1)+1, 'MarkerFaceColor','r'); hold on; 
                    plot(xm,ym,'ko','MarkerSize',sqrt(jj+1)+1); hold on;              
                    axis equal;
                    title('ori pair walk');        
                subplot(2,2,3);
                    plot(0*oridist_table(:,1)+jj,oridist_table(:,2),'bo','MarkerSize',5, 'MarkerFaceColor','b'); hold on;
                    plot(0*oridist_table(:,1)+jj,oridist_table(:,3),'ro','MarkerSize',5, 'MarkerFaceColor','r'); hold on;
                    plot(0*oridist_table(:,1)+jj,oridist_table(:,4),'ro','MarkerSize',7, 'MarkerFaceColor','w'); hold on; 
                    plot(0*oridist_table(:,1)+jj,oridist_table(:,5),'ks','MarkerSize',7, 'MarkerFaceColor','w'); hold on;     
                    %legend('ter-cellcom'   ,      'ori-nearori'     ,    'oricom-ter'  ,      'oricom-cellcom',...
                    %    'Location','SouthOutside');
                    xlabel('frame index');
                    ylabel('d/minor axis, r.u.'); 
                    
            function Plot_replisome_panels(MatFilePath,CellFrameName,figureno);
                %plot ter-ori relations
                %repli_pair: [movieframe xo xo2 xp yo yo2 yp xm ym d_oo
                figure(1);
                CellInfoPath=strcat(MatFilePath,'ResultsOfCell',CellFrameName);
                load(strcat(CellInfoPath,'_SpotPairs.mat'),'replipair_table');
                load(strcat(CellInfoPath,'_Clusters.mat'),'Clusters','chro_pic','celledge_pic');
                load(strcat(CellInfoPath,'_Cellshape.mat'),'Cell');    
                load(strcat(CellInfoPath,'_Spots.mat'), 'All_labels');
                frames=replipair_table(:,1);
                distance=replipair_table(:,10);
                plot(frames,distance,'o'); hold on;
                