function A080_WF_BuildDistanceDensityMaps(batchrunindex)
%JWJK_A:-------------------------------------------------------------------
%
%Building of Dynamic 'Hi-C' analogons
%
%Summary: This function build 'Hi-C' style 2D-plots using the circular mapping 
%of the expanded chromatin 
%
%Approach: We use the relation between spatial distance and genomic content 
%along the ridge of the donut-shaped chromatin pattern (the '1D density curve').
%as with an Hi-C map, the axes are genomic contact. as a measure for
%contact, or mutual vicinity of genomic material, we use the negative of the
%spatial distance as parameter (i.e., shorted distance means highest 'contact' signal)
%and plot that in 2D.
%
%Input: data in .mat files stored in former analysis steps.
%
%Output: Data is saved as images. 
%
%:JWJK_A-------------------------------------------------------------------

close all;

JustMakeControlKymo=0;
ExampleIndex=1;
ControlKymo=[];

if nargin<1,batchrunindex=11.1; end
initval=A000__WF_Get_JacobPathsandExperiments(batchrunindex);
PadPerc=initval.Padcurves;


FewCells=0;  %Demo &snapshot purposes
if FewCells
    initval.Cell_Labels=[{'1000002_t09'}; {'1000002_t09'}];  %11.1
%      initval.Cell_Labels=[{'200020'}; {'200022'}];           %62
%      initval.Cell_Labels=[{'100578'}; {'100578'}];           %9
end

allframes=length(initval.Cell_Labels);
imoutdir=strcat(initval.resultpath,'A080_ContactMaps',initval.DirSep);
if ((~JustMakeControlKymo)&(~FewCells));
    imoutdir=strcat(initval.resultpath,'A080_ContactMaps',initval.DirSep);
    
%     if isdir(imoutdir), rmdir(imoutdir,'s');  end
%     mkdir(imoutdir);
%     mkdir(strcat(imoutdir,initval.DirSep,'Tiffs'));
%     mkdir(strcat(imoutdir,initval.DirSep,'Plots'));
    
end
cnt=0;
for jj=1:allframes
    disp(strcat(num2str(allframes-jj+1),'Contact maps to go'));
    cellno=char(initval.Cell_Labels{jj});
    CellName=strcat('ResultsOfCell',cellno,'.mat'); 
    MatFilePath=strcat(initval.resultpath,'ResultsPerCellMatlab',initval.DirSep);
    load(strcat(MatFilePath,CellName));
    if 0
       MX=Chromosome.CartesianContourMax_X;
       MY=Chromosome.CartesianContourMax_Y;
       pcolor(Chromosome.picture); colormap bone; shading flat; hold on;
       plot(MX,MY,'y-','LineWidth',2); hold on;
       plot(Cfp.spotX,Cfp.spotY, 'bo','MarkerSize',10, 'MarkerFaceColor','b'); hold on;
       plot(Rfp.spotX,Rfp.spotY, 'ro','MarkerSize',10,'MarkerFaceColor','r'); hold on;
       [~]=ginput(1);
    end
    
    if GeneralCellProps.Okayproduct
        cnt=cnt+1;
        
        %Swap the equidistance from  distance axis to bpair axis 
        DistAx=Aligned.Dist.NormAxis;
        NormCumDensity_vs_Dist_ori=Aligned.Dist.NormCumDensityMarkercorrected;
        
        %First, build local density curve
        Density=diff(NormCumDensity_vs_Dist_ori);
        Density=[Density(end) Density];  
        
        
        
        RunDirection='CW';
        switch RunDirection
            case 'CW' %choose: density in CW direction     
             NormCumDensity_vs_Dist=cumsum(Density);
            case 'CCW' %choose: density in CCW direction     
            NormCumDensity_vs_Dist=cumsum(fliplr(Density));
        end
        
        dum=1;
        [DistAx,NormCumDensity_vs_Dist]=MakeMonotonicallyRising(DistAx,NormCumDensity_vs_Dist);            
        
        Bp_Normax=1:100;
        
        %re-interpolate vfor equidistant points
        if (length(NormCumDensity_vs_Dist)>2)
        Distance_vs_BP=interp1(NormCumDensity_vs_Dist,DistAx,Bp_Normax,'Linear',0);  
        %terpos=Aligned.Orig.MarkerLabel.MarkerBPPercExpected;   
        
        
        [Bp_Normax,Distance_vs_BP]=Pad_it(Bp_Normax,Distance_vs_BP,25);
        
        
 
        LB=length(Distance_vs_BP);
      
       %build and save the heat map per image
       if ~JustMakeControlKymo
           DistanceHeatMap=zeros(LB,LB); 
            for ii=1:LB
                for jj=1:LB
                    distance=abs(Distance_vs_BP(ii)-Distance_vs_BP(jj));
                    DistanceHeatMap(ii,jj)=distance;
                end
            end 
            
            %build an averaged heatmap
            if cnt==1
                AverageDistanceHeatMap=DistanceHeatMap;
            else
                AverageDistanceHeatMap=AverageDistanceHeatMap... 
                    +DistanceHeatMap;
            end


            %set the cutoff
            if 0
            MX=max(DistanceHeatMap(:));
            MN=min(DistanceHeatMap(:));
            tres=0.2*MX;
            sel=find(DistanceHeatMap>tres);
            DistanceHeatMap(sel)=tres;
            end
            pcolor((-DistanceHeatMap)); shading flat; colormap hot; hold on;

            [rr,cc]=size(DistanceHeatMap);


            %plot the result
            plot(1:LB,(1:LB)*0+LB/2,'k-');
            plot((1:LB)*0+LB/2,1:LB,'k-');
            xlabel('Approximate genomic Distance, %');
            ylabel('Approximate genomic Distance, %');
            axis square
            axis ij

            OutName=strcat('DistanceMap',CellName(1:end-4));

            %build a tiff picture
            imout=uint16(10*(100-DistanceHeatMap-1));
            imwrite(imout,strcat(imoutdir,'Tiffs',initval.DirSep,OutName,'.tif'),'tif');
            if 1
            saveas(gcf,strcat(imoutdir,'Plots',initval.DirSep, OutName,'.jpg')); 
            pause(0.01);
            end
            hold off;
       else %just build a control kymograph
           crv=diff(Distance_vs_BP);
           ControlKymo=[ControlKymo; crv];
           if cnt==ExampleIndex, plot_Distance_vs_BP=Distance_vs_BP; end
       end
        end
    end
end
if JustMakeControlKymo
    figure;
    subplot(2,2,1);
    plot(Bp_Normax,plot_Distance_vs_BP,'k-','Linewidth', 2);
    title('Distance% vs. BP%')
    xlabel('basepair %')
    ylabel('distance %');
    text(10,95,CellName(10:end-4));
    subplot(2,2,3);
    plot(Bp_Normax(2:end),diff(plot_Distance_vs_BP),'k-','Linewidth', 2);
    title('Gaps: Derivative Distance% vs. BP%')
    xlabel('basepair %')
    ylabel('d(distance)/d(BP)');
    
    subplot(1,2,2);
    ControlKymo=(max(ControlKymo(:))-ControlKymo);
    pcolor(ControlKymo); shading flat, colormap jet
    xlabel('basepair %')
    ylabel('movie frames sorted by cell');
    title(strcat('Gap map of ',initval.expi));
    saveas(gcf,strcat(initval.resultpath,initval.DirSep,initval.expi,'A080_GapMap.jpg'),'jpg');
    saveas(gcf,strcat(initval.resultpath,initval.DirSep,initval.expi,'A080_GapMap'));
    pause(0.01);
else



%% process, plot and save average
    AverageDistanceHeatMap=AverageDistanceHeatMap/cnt;
    
    if 0
    MX=max(AverageDistanceHeatMap(:));
    MN=min(AverageDistanceHeatMap(:));
    tres=0.2*MX;
    sel=find(AverageDistanceHeatMap>tres);
    AverageDistanceHeatMap(sel)=tres;
    end   
 figure;   
    pcolor((-AverageDistanceHeatMap)); shading flat; colormap hot; hold on;
    plot(1:LB,(1:LB)*0+LB/2,'k-');
    plot((1:LB)*0+LB/2,1:LB,'k-');
    xlabel('Approximate genomic Distance, %');
    ylabel('Approximate genomic Distance, %');
    axis square
    axis ij

    
    %build a tiff picture
    imout=uint16(10*(100-AverageDistanceHeatMap-1));
    imwrite(imout,strcat(imoutdir,'AverageDistanceMapOf',initval.expi,'.tif'),'tif');
    saveas(gcf,strcat(imoutdir,'AverageDistanceMapOf',initval.expi,'.jpg')); 
    pause(0.1);
    hold off;
end     
    
    function [bbax_out,bb_out]=MakeMonotonicallyRising(bbax_in,bb_in)
        %function makes flat or negative parts -slightly- rising
        bb_out=0*bb_in;
        Lbb=length(bb_in);
        cnt=1;
        bb_out(1)=bb_in(1); bb_out(end)=bb_in(end);
        for ii=1:Lbb-1
            bb_out(ii+1)=nanmax([nanmax(bb_in(1:ii)) bb_in(ii+1)]);
        end
        [bb_out,idx]=unique(bb_out);
        bbax_out=bbax_in(idx);
        %plot(bbax_in,bb_in,'b-'); hold on;
        %plot(bbax_out,bb_out,'ro-');
        sel=~isnan(bb_out); bb_out=bb_out(sel); bbax_out=bbax_out(sel);
        dum=1;
        
        function [Axzout, Crvout]=Pad_it(Axin,Crvin,Perc);
            %This function pads Dist-density curves 
            LL=length(Axin); CL=ceil(LL*Perc/100);
            PreAx=Axin(1:CL); PostAx=Axin(end-CL:end);
            PreCrv=Crvin(1:CL); PostCrv=Crvin(end-CL:end);
            
            Pre2PostAx=PreAx-PreAx(1)+PostAx(end)+PreAx(2)-PreAx(1);
            Post2PreAx=PostAx-PostAx(end)+PreAx(1)-(PostAx(end)-PostAx(end-1));
            
            Pre2PostCrv=PreCrv-PreCrv(1)+PostCrv(end)+PreCrv(2)-PreCrv(1);
            Post2PreCrv=PostCrv-PostCrv(end)+PreCrv(1)-(PostCrv(end)-PostCrv(end-1));
            
            Axzout=[Post2PreAx Axin Pre2PostAx];
            Crvout=[Post2PreCrv Crvin Pre2PostCrv];
            dum=1;
            

