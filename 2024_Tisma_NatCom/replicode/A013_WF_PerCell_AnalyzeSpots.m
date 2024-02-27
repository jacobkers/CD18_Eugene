function A013_WF_PerCell_AnalyzeSpots(initval)
%JWJK_A:-------------------------------------------------------------------
%Description: analyze local lables (typically ori and ter) in their
%respective channels

%input: .mat database from A010 analysis.

%output:.mat database is expanded, cell images are saved; some tables. All
%with 'A013' label.

%Reference: CD lab, project Sandro, written by Jacob Kers 2018-20
%:JWJK_A-------------------------------------------------------------------

% c1_cell_001t1xy1.mat contains cellc1 stack: cfp
% c3_cell_001t1xy1.mat contains cellc3 stack: rfp
% cellc3

close all;
%runtime options
actions.showandsavepics=1;
if nargin<1
    usr='Jacob', batchrunindex=0;
    initval=A000_Repli_Init(batchrunindex,usr);
end

disp(initval.expi);
CellImagePath=strcat(initval.pth_repli,'CellImages_All',initval.DirSep);
if isdir(CellImagePath), 
    rmdir(CellImagePath,'s');  
end
mkdir(CellImagePath);
MatFilePath=strcat(initval.pth_repli,'ResultsPerCellMatlab',initval.DirSep);
LC=length(initval.Cell_Labels);



for ii=1:LC
    cellno=char(initval.Cell_Labels{ii});    
    CellName=strcat('ResultsOfCell',cellno); 
    disp(strcat('Program:A13_experiment:',initval.expi,':',num2str(LC-ii), 'to go-',CellName,'Label localization..', num2str(LC-ii+1), 'cells to go'));
    %cellno=num2str(initval.cell_labelnos(ii));    
    NumCellLabel=BuildNumericCellLabel(cellno);
    CellSpecs=[ii NumCellLabel];   
    
    channel_stack=get_channel_stack_from_cropcode(cellno,initval);
    [ cell_pic, cellmask, chro_pic, Rfp_pic, Cfp_pic, ]=get_channel_ID(channel_stack, initval);

    %bit of refinement
    chro_pic=(chro_pic-min(chro_pic(:))).*cellmask;
    edge_pic = bwmorph(cellmask,'remove'); 
    
    CellshapeName=strcat(MatFilePath,strcat('ResultsOfCell',cellno,'_Cellshape.mat'));
    load(CellshapeName, 'Cell');
%% 
   %------------------------------------------------------------
   %3)The respective ter and ori images ('cfp' and 'rfp') are loaded.  
   % For each, a single spot center and corresponding spot properties 
   % (position, intensity) is tracked following 
   % Llorente-Garcia / Reyes et al.
   if initval.analyze_ori_ter
        Cfp=Work_spotimage_CfpTer(Cfp_pic,cellmask,Cell, initval);
        All_labels.Cfp=Cfp;
        All_labels.Cfp_pic=Cfp_pic;   
    end
    if initval.analyze_replisome | initval.analyze_ori_ter
        Rfp=Work_spotimage_RfpOri(Rfp_pic, cellmask,Cell, initval);
        Rfp.frameindex=0*Rfp.spotX+ii;
        All_labels.Rfp=Rfp;
        All_labels.Rfp_pic=Rfp_pic;  
    end
    
    
    %save per Cell  
    initval.ResultName=strcat(MatFilePath,strcat('ResultsOfCell',cellno,'_Spots.mat'));        
    save(initval.ResultName, 'All_labels');

    %collect
    if initval.analyze_replisome
        replisummary.N_all(ii)=length(Rfp.spotX)*Rfp.OK_nicepic(1);
        replisummary.N_strong(ii)=sum(Rfp.OK_strong)*Rfp.OK_nicepic(1);
        replisummary.N_pair(ii)=sum(Rfp.OK_pair)*Rfp.OK_nicepic(1);
        replisummary.N_pair(ii)=sum(Rfp.OK_pair)*Rfp.OK_nicepic(1);
        replisummary.N_goodpic(ii)=Rfp.OK_nicepic(1);
    end
    
    
    
    if actions.showandsavepics
         if initval.showplots
            set(figure(1), 'visible','on');
         else
            set(figure(1), 'visible','off');
         end      
           subplot(2,2,1);
               edge_sel=find(edge_pic); plotpic=chro_pic; plotpic(edge_sel)=max(plotpic(:)); 
               pcolor(plotpic); shading flat, colormap bone; hold on; axis equal;            
               if initval.analyze_replisome||initval.analyze_ori_ter
                    %weakspots
                    wk=find(Rfp.OK_strong==0);
                    if length(wk)>0,
                    plot(Rfp.spotX(wk),Rfp.spotY(wk), 'yx', 'MarkerSize', 6); 
                    end
                    %strongspots
                    str=find(Rfp.OK_strong==1);
                    plot(Rfp.spotX(str),Rfp.spotY(str), 'wo', 'MarkerSize', 8); 
                    %build pairs
                    sel=find(Rfp.OK_pair==1);
                    LP=length(sel); pairX=zeros(LP,2); pairY=zeros(LP,2);
                    for pp=1:LP
                        id1=sel(pp); 
                        id2=Rfp.OK_pair_idx(id1);
                        pairX(pp,:)=[Rfp.spotX(id1) Rfp.spotX(id2)];
                        pairY(pp,:)=[Rfp.spotY(id1) Rfp.spotY(id2)];           
                    end       
                    plot(pairX',pairY','wo-', 'MarkerSize',4,'MarkerFaceColor','w');
                    %text(5,5,['#',num2str(Cell.N_est)],'Color','white');                  
               end            
               end
                if initval.analyze_ori_ter
                         plot(Cfp.spotX,Cfp.spotY, 'bo','MarkerSize',10); hold on;                        
                end 
                title(strcat('Cell', num2str(cellno,'% 3.0f')));
                pause(0.01);
            
           if initval.analyze_replisome||initval.analyze_ori_ter
           subplot(2,2,2);
                    Rfp_plotpic=Rfp_pic;
                    Rfp_plotpic(edge_sel)=max(Rfp_plotpic(:)); 
                    pcolor(Rfp_plotpic); shading flat, colormap bone; axis equal; axis off; hold on;   
                                          %weakspots
                        wk=find(Rfp.OK_strong==0);
                        if length(wk)>0,
                        plot(Rfp.spotX(wk),Rfp.spotY(wk), 'yx', 'MarkerSize', 6); 
                        end
                        %strongspots
                        str=find(Rfp.OK_strong==1);
                        plot(Rfp.spotX(str),Rfp.spotY(str), 'wo', 'MarkerSize', 8); 
                        %build pairs
                        sel=find(Rfp.OK_pair==1);
                        LP=length(sel); pairX=zeros(LP,2); pairY=zeros(LP,2);
                        for pp=1:LP
                            id1=sel(pp); 
                            id2=Rfp.OK_pair_idx(id1);
                            pairX(pp,:)=[Rfp.spotX(id1) Rfp.spotX(id2)];
                            pairY(pp,:)=[Rfp.spotY(id1) Rfp.spotY(id2)];           
                        end       
                        plot(pairX',pairY','wo-', 'MarkerSize',4,'MarkerFaceColor','w');
                    title('rfp');
           end
                
           if initval.analyze_ori_ter
                subplot(2,2,3);
                    Cfp_plotpic=Cfp_pic;
                    Cfp_plotpic(edge_sel)=max(Cfp_plotpic(:)); 
                    pcolor(Cfp_plotpic); shading flat, colormap bone; hold on;
                    plot(Cfp.spotX,Cfp.spotY, 'bo','MarkerSize',10); hold on;        
                     title('cfp');
           end    
            pause(0.01);           
            saveas(gcf,strcat(initval.pth_repli,'CellImages_All',initval.DirSep,'Cell', num2str(cellno,'% 3.0f'),'.jpg')); 
            pause(0.01); 
            %[~]=ginput(1); 
            close(gcf);
end

if initval.analyze_replisome %summary
    RepliSummaryLegend=([{'1.initialcount'},{'2.strengthOK'},{'3.paired'},{'4.#ofpicsused'}]);
    RepliSummaryVals=[nansum(replisummary.N_all)  nansum(replisummary.N_strong)...
                      nansum(replisummary.N_pair) nansum(replisummary.N_goodpic)];

    LS=length(RepliSummaryVals);
    axz=1:LS;
    plot(axz,0*axz,'k-'); hold on;
    for ii=1:LS
        bar(axz(ii),RepliSummaryVals(ii)); hold on;    
    end
    dum=1;
    legend(RepliSummaryLegend,'Location','NorthOutside');
    title(strcat('SpotScreening Report of: ',initval.expi));
    saveas(gcf,strcat(initval.pth_repli,initval.DirSep,initval.expi,'_A013_SpotScreeningReport.jpg'),'jpg');
end

