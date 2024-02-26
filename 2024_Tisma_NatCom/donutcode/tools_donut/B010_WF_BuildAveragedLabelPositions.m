function B010_WF_BuildAveragedLabelPositions(batchrunindex)
%JWJK_B:-------------------------------------------------------------------
%Title: Get and store averaged averaged label position
%
%Summary: for movies, define an averaged albel position for a more steady
%monitoring of cluster dynamics
%
%Input: data in .mat files stored in former analysis steps.
%
%Output: 'Cfp_Av, Rfp_Av structures
%
%:JWJK_B-------------------------------------------------------------------
sho=0;
close all;
CellNames=cell(0);
if nargin<1,
    batchrunindex=11;
    sho=1; 
end
initval=A000__WF_Get_JacobPathsandExperiments(batchrunindex);
disp(initval.expi);    


%% After doing all cells, check if database was a movie. If so, obtain averaged cfp and rfp positions and save them
%First, make a list of cells via the first frame
sel=strfind(initval.Cell_Labels,'_t01');

for ii=1:length(sel);
    if ~isempty(sel{ii})
        CellNames=[CellNames initval.Cell_Labels{ii}];
    end
end
NCells=length(CellNames);
for cc1=1:NCells  %for all cells:
    disp(strcat('Getting Average Label Posses..',num2str(NCells-cc1+1), 'cells to go'));
    CellName=char(CellNames(cc1));
    CellName=CellName(1:end-4);
    CellFrames=cell(0);    
    sel=strfind(initval.Cell_Labels,CellName);
    for ii=1:length(sel);
        if ~isempty(sel{ii})
            CellFrames=[CellFrames initval.Cell_Labels{ii}];
        end
    end
    
    
    
    Frs1=length(CellFrames);   
    RadPosses_rfp=zeros(Frs1,1);
    RadPosses_cfp=zeros(Frs1,1);
    
    for fi=1:Frs1  %for all frames of a cell; collect spot properties
        cellno1=char(CellFrames{fi});
        FrameTime=str2num(cellno1(end-1:end));
        CellName1=strcat('ResultsOfCell',cellno1,'.mat'); 
        MatFilePath=strcat(initval.resultpath,'ResultsPerCellMatlab',initval.DirSep);
        load(strcat(MatFilePath,CellName1),'Cfp_pic','Rfp_pic','Cfp','Rfp','Chromosome');
        
                
        %load the stuff
        Cfp_pic=double(Cfp_pic);
        Rfp_pic=double(Rfp_pic);
        
        extraedge=1; 
        switch initval.cellmaskname  %SuperUglyHack            
        case 'cellmask'; 
        load(strcat(initval.maindatapath,'ma_',cellno1,'.mat'), initval.cellmaskname);
        cellmask=1.0*bwmorph(cellmask,'dilate',extraedge);
        case 'cellma'
          load(strcat(initval.maindatapath,'ma_',cellno1,'.mat'), 'cellma');
          cellmask=1.0*bwmorph(cellma,'dilate',extraedge);
        end
        maskarea=find(cellmask==0); 
        cellarea=find(cellmask>0); 

        Cfp_pic(maskarea)=min(Cfp_pic(cellarea));
        Rfp_pic(maskarea)=min(Rfp_pic(cellarea));
        
         if fi==1
            Cfp_pic_av=Cfp_pic;
            Rfp_pic_av=Rfp_pic;
         else
             Cfp_pic_av=Add2AvPic(Cfp_pic_av,Cfp_pic);
             Rfp_pic_av=Add2AvPic(Rfp_pic_av,Rfp_pic);
         end       
         RadPosses_cfp(fi)=Cfp.spotPosAngle;
         RadPosses_rfp(fi)=Rfp.spotPosAngle;
         
                 
        if sho
            figure(1);
            subplot(2,2,1); pcolor(Cfp_pic); shading flat;
            title('Cfp-snapshot');
            subplot(2,2,3); pcolor(Rfp_pic); shading flat;   
            title('Rfp-snapshot');
            subplot(2,2,2); pcolor(Cfp_pic_av); shading flat; 
            title('Cfp-average');
            subplot(2,2,4); pcolor(Rfp_pic_av); shading flat; 
            title('Rfp-average');
            pause(0.1);
        end       
    end    
    %some cleaning;    
    Cfp_pic_av=Cfp_pic_av/Frs1;
    maskarea=find(isnan(Cfp_pic_av)); 
    cellarea=find(~isnan(Cfp_pic_av)); 
    Cfp_pic_av(maskarea)=min(Cfp_pic_av(cellarea));
   
    %some cleaning;  
    Rfp_pic_av=Rfp_pic_av/Frs1;
    maskarea=find(isnan(Rfp_pic_av)); 
    cellarea=find(~isnan(Rfp_pic_av)); 
    Rfp_pic_av(maskarea)=min(Rfp_pic_av(cellarea));
        
    %get spot properties
    Cfp_AvProps=Get_SpotProps(Cfp_pic_av,Chromosome,initval);
    Rfp_AvProps=Get_SpotProps(Rfp_pic_av,Chromosome,initval);
    
    
    
    if sho
        figure(3);
        plot(RadPosses_cfp,'bo-'); hold on;        
        plot(0*RadPosses_cfp+median(RadPosses_cfp),'bo-'); hold on;
        plot(0*RadPosses_cfp+Cfp_AvProps.spotPosAngle,'b-'); hold on;   
        plot(RadPosses_rfp,'ro-'); hold on;
        plot(0*RadPosses_rfp+median(RadPosses_rfp),'ro-'); hold on;
        plot(0*RadPosses_rfp+Rfp_AvProps.spotPosAngle,'r-'); hold on;                
        
        title('Radial position')
        legend('cfp per pic','cfp-median','cfp from averaged pic',...
               'rfp per pic','rfp-median','rfp from averaged pic');
        xlabel('frames');
        ylabel('rads');   
        [~]=ginput(1);
        close(gcf);
    end
            
    for fi=1:Frs1  %for all frames of a cell
        %re_run the cell frames to save the avareaged result to every
        %single cell frame
        cellno1=char(CellFrames{fi});
        FrameTime=str2num(cellno1(end-1:end));
        CellName1=strcat('ResultsOfCell',cellno1,'.mat'); 
        MatFilePath=strcat(initval.resultpath,'ResultsPerCellMatlab',initval.DirSep);
        initval.ResultName=strcat(MatFilePath,strcat('ResultsOfCell',cellno1,'.mat'));
        save(initval.ResultName,'Cfp_AvProps','Rfp_AvProps','-append');
        %dum=1;     
    end
end
dum=1;

function pic_av_out=Add2AvPic(pic_av_in,pic);      
        pic=pic-nanmedian(pic(:));
        pic_av_in=pic_av_in-nanmedian(pic_av_in(:));
        
       
        [rr0,cc0]=size(pic);
        [rr1,cc1]=size(pic_av_in);
        rr_cr=min([rr0 rr1]);
        cc_cr=min([cc0 cc1]);
        %sync sizes
        pic=pic(1:rr_cr,1:cc_cr);
        pic_av_in=pic_av_in(1:rr_cr,1:cc_cr);
        
        %add up
        pic_av_out=pic_av_in+pic;
%         figure(2); plot(sum(pic_av_in)); hold on;           
%         dum=1;
        
        

