function A062_MukB_Dna_overlap(initval)
%JWJK_A:-----------------------------------------------------------------
%Description: process colocalization data 

%input: 
    %1) .mat databases from A010/013/060 replicode analysis.
    %2) .mat database from a specified directory containing donut analysis
    %results

%output: directory labeled A062 with images, table, database

%Reference: CD lab, project Sandro, written by Jacob Kers 2018-20
%:JWJK_A-------------------------------------------------------------------

close all;
save_jpg=1;
axislimit=6;

if nargin<1
    usr='Jacob', 
    batchrunindex=24;
    override_paths=0;
    initval=A000_Repli_Init(batchrunindex,usr);
end


Psf_meas=initval.Psf_est;
allframes=length(initval.Cell_Labels);


UsedCellNames=[]; ExcludedCellNames=[]; used_jj=0;
%load tresholded images
load([initval.pth_repli, '\',initval.expi, '_A060_AllCellsResults.mat'], ...
                     'Allcells_Mukbef_Chro_images',...
                     'Allcells_Ori_Chro_images',...
                     'Allcells_Ori_Ter_images',...
                     'UsedCellNames',... 
                     'ExcludedCellNames');
                 
imoutdir=strcat(initval.pth_repli,'A062_tresholds',initval.DirSep);  
if isdir(imoutdir)
    rmdir(imoutdir,'s'); 
end
mkdir(imoutdir);   

celloutnames=[];
expoutnames=[];
symmetries=[];

% run the cells
for jj=1:allframes     
    cellno=char(initval.Cell_Labels{jj}); 
    CellName=strcat('ResultsOfCell',cellno);
    disp(strcat('Program:A62_experiment:',initval.expi,':',CellName,'MukB_DNA_overlap_analysis..', num2str(allframes-jj+1), 'cells to go'));   
    notexcluded=CheckUserExclusion(cellno, initval);
    if notexcluded
        celloutnames=[celloutnames; {CellName}];
        expoutnames=[expoutnames; {initval.expi}];
        %load per Cell  
        ThisCellLoadName=strcat(initval.pth_repli,'\ResultsPerCellMatlab\',CellName,'_Spots.mat');        
        load(ThisCellLoadName, 'All_labels'); 
        ThisCellLoadName2=strcat(initval.pth_repli,'\ResultsPerCellMatlab\',CellName,'_Cellshape.mat'); 
        load(ThisCellLoadName2); 
    
    for ii=1:initval.numberofchannels
        [refpic,cellmask]=Get_channel_pic(initval,cellno,2); %chromosome
        [pic,~]=Get_channel_pic(initval,cellno,ii);
        if ii==1
            [rr,cc]=size(pic);
            picstack=zeros(rr,cc,initval.numberofchannels);
        end
        picstack(:,:,ii)=pic;
    end
    
    % collect data we need
    %note some swapping coordinates etc. to get it right
    im_mukb=picstack(:,:,5)';
    im_dna=picstack(:,:,2)';
    im_ori=picstack(:,:,3)';
    im_ter=picstack(:,:,4)';
    
    %find properindex
    used_jj=find(strcmp(cellno, UsedCellNames)==1);
    
    im_tres_mukb=Allcells_Mukbef_Chro_images(used_jj).masked_image_A;
    im_tres_dna=Allcells_Mukbef_Chro_images(used_jj).masked_image_B;
    im_tres_ori=Allcells_Ori_Chro_images(used_jj).masked_image_A;
    im_tres_ter=Allcells_Ori_Ter_images(used_jj).masked_image_B;
    %ori pos:
    xo=All_labels.Rfp.spotY(1);
    yo=All_labels.Rfp.spotX(1);
    %donut com pos:
    xc=Cell.Centroid(2);
    yc=Cell.Centroid(1);
    %coordinates
    [rr,cc]=size(im_dna);
    [XX,YY]=meshgrid(1:cc,1:rr);
    
    %add an extra circular mask around ori
    RR=((XX-xo).^2+(YY-yo).^2).^0.5;
    ori_mask=1.0*RR<initval.orimaskradius;
    ori_mask_edge= bwmorph(ori_mask,'remove');
    edge_xx=XX(ori_mask_edge>0);
    edge_yy=YY(ori_mask_edge>0);
    
    % Get symmetry of pattern relative to coordinates
    %use(im_dna, xc,yc,xo,yo);
    %out: symmetry,
    %first, get rotation 
    
    alpha=90-180/pi*atan2(yo-yc,xo-xc);
    [xor,yor]=Rotate_Points(xc,yc,xo,yo,alpha);
    axislength=((xor-xc)^2+(yor-yc)^2)^0.5;
    
    %rotate dna
    ix_dna=find(im_tres_dna>0 & ori_mask>0);
    area_dna=length(ix_dna);
    if area_dna>0
        %rotate dna-------------------
        %2D coordinates
        xx_dna=XX(ix_dna);
        yy_dna=YY(ix_dna);    
        [xx_dna_r,yy_dna_r]=Rotate_Points(xc,yc,xx_dna,yy_dna,alpha);
        %split in left and right of now-vertical axis
        sel_dna_left=find(xx_dna_r<=xc);
        sel_dna_right=find(xx_dna_r>xc);          
        %symmetry area-only
        area_dna_left=length(sel_dna_left);
        area_dna_right=length(sel_dna_right);
        symmetry_dna_by_area=(area_dna_left-area_dna_right)/area_dna;
        %image indices
        ix_dna_left=ix_dna(sel_dna_left);
        ix_dna_right=ix_dna(sel_dna_right);
        %symmetry by intensity
        intensity_dna=sum(im_dna(ix_dna));
        intensity_dna_left=sum(im_dna(ix_dna_left));
        intensity_dna_right=sum(im_dna(ix_dna_right));
        symmetry_dna_by_intensity=(intensity_dna_left-intensity_dna_right)/intensity_dna;
        
        %same, mukb--------------------
        ix_mukb=find(im_tres_mukb>0 & ori_mask>0);
        area_mukb=length(ix_mukb);
        %2D coordinates
        xx_mukb=XX(ix_mukb);
        yy_mukb=YY(ix_mukb);    
        [xx_mukb_r,yy_mukb_r]=Rotate_Points(xc,yc,xx_mukb,yy_mukb,alpha);
        %split in left and right of now-vertical axis
        sel_mukb_left=find(xx_mukb_r<=xc);
        sel_mukb_right=find(xx_mukb_r>xc);
        %symmetry area-only
        area_mukb_left=length(sel_mukb_left);
        area_mukb_right=length(sel_mukb_right);
        symmetry_mukb_by_area=(area_mukb_left-area_mukb_right)/area_mukb;
        %image indices
        ix_mukb_left=ix_mukb(sel_mukb_left);
        ix_mukb_right=ix_mukb(sel_mukb_right);
        %symmetry by intensity
        intensity_mukb=sum(im_mukb(ix_mukb));
        intensity_mukb_left=sum(im_mukb(ix_mukb_left));
        intensity_mukb_right=sum(im_mukb(ix_mukb_right));
        symmetry_mukb_by_intensity=(intensity_mukb_left-intensity_mukb_right)/intensity_mukb;
                        
        %collect
        symmetries=[symmetries; [symmetry_dna_by_area, symmetry_mukb_by_area, symmetry_dna_by_intensity, symmetry_mukb_by_intensity, axislength, 1.0*(axislength>axislimit)] 1.0*(area_dna>area_mukb)];    
        ok_treshold=(1.0*(axislength>initval.axislimit) && 1.0*(area_dna>area_mukb));
    end
    
    % plot menu
    if 1 
         set(figure(1), 'visible','off');
        %check xy-pix etc
        subplot(2,4,1); pcolor(im_dna); shading flat; colormap jet; %dna 
        title('DNA');
        axis equal;
        subplot(2,4,2); pcolor(im_ori); shading flat; colormap jet; %ori
        title('Ori'); axis equal;
        subplot(2,4,3); pcolor(im_ter); shading flat; colormap jet; %ter
        title('Ter'); axis equal;
        subplot(2,4,4); pcolor(im_mukb); shading flat; colormap jet; %mukb
        title('MukB'); axis equal;      
        subplot(2,4,5);
         pcolor(im_tres_dna); shading flat; colormap jet; hold on;
         plot(xc,yc,'wo', 'MarkerSize', 4,'MarkerFaceColor','w');
         plot(xo,yo,'ro', 'MarkerSize', 4, 'MarkerFaceColor','r');
          plot([xo xc],[yo yc],'w-');
          plot(edge_xx,edge_yy, 'wo','MarkerSize', 1);
          axis equal;
         subplot(2,4,6);
         pcolor(im_tres_ori); shading flat; colormap jet; hold on;
         axis equal;
         subplot(2,4,7);
         pcolor(im_tres_ter); shading flat; colormap jet; hold on;
         axis equal;         
       subplot(2,4,8);
         pcolor(im_tres_mukb); shading flat; colormap jet; hold on;
          plot(xc,yc,'wo', 'MarkerSize', 4,'MarkerFaceColor','w');
         plot(xo,yo,'ro', 'MarkerSize', 4, 'MarkerFaceColor','r');
          plot([xo xc],[yo yc],'w-');
          plot(edge_xx,edge_yy, 'wo','MarkerSize', 1);
         axis equal;
         
         imoutdir=strcat(initval.pth_repli,'A062_tresholds\');
         if ok_treshold
            outname=['accepted_',CellName, '_A062_tresholds_overview'];  
         else
            outname=['rejected_',CellName, '_A062_tresholds_overview'];
         end
         if save_jpg            
            saveas(gcf,[imoutdir, outname, '_channels.jpg' ]); 
         end

        set(figure(2), 'visible','off')
        subplot(2,2,1); 
            plot(xx_dna,yy_dna,'ko', 'MarkerFaceColor', 'k'); hold on;
            %plot([xo xc],[yo yc],'k-');
            title('DNA');
            axis equal;
        subplot(2,2,3); 
            plot(xx_dna_r(sel_dna_left),yy_dna_r(sel_dna_left),  'ro', 'MarkerFaceColor', 'r'); hold on
            plot(xx_dna_r(sel_dna_right),yy_dna_r(sel_dna_right),'bo','MarkerFaceColor', 'b'); hold on
            %plot([xor xc],[yor yc],'k-');
            %text(xc+1,yc,num2str(symmetry_dna));
            %text(xc+1,yc+5,num2str(axislength));
            axis equal;
        subplot(2,2,2); 
            plot(xx_mukb,yy_mukb,'ko','MarkerFaceColor', 'k'); hold on;
            %plot([xo xc],[yo yc],'k-');
            title('MukB');
            axis equal;
        subplot(2,2,4); 
            plot(xx_mukb_r(sel_mukb_left),yy_mukb_r(sel_mukb_left),'ro', 'MarkerFaceColor', 'r'); hold on
            plot(xx_mukb_r(sel_mukb_right),yy_mukb_r(sel_mukb_right),'bo','MarkerFaceColor', 'b'); hold on
            %plot([xor xc],[yor yc],'k-');
            %text(xc+1,yc,num2str(symmetry_mukb));
            axis equal;
            
         if save_jpg
            saveas(gcf,[imoutdir, outname, '_rotations.jpg' ]);
%           else
%             [~]=ginput(1); 
         end
         close all;       
    end

    end
end

%process good ones
okcells=find((symmetries(:,6)==1) & (symmetries(:,7)==1));
N_ok=length(okcells);
%by area
symm_dna_area=abs(symmetries(okcells,1));
symm_mukbef_area=abs(symmetries(okcells,2));
symm_dna_md_area=nanmedian(symm_dna_area);
symm_dna_mn_area=nanmean(symm_dna_area);
symm_mukbef_md_area=nanmedian(symm_mukbef_area);
symm_mukbef_mn_area=nanmean(symm_mukbef_area);
%by intensity
symm_dna_intensity=abs(symmetries(okcells,3));
symm_mukbef_intensity=abs(symmetries(okcells,4));
symm_dna_md_intensity=nanmedian(symm_dna_intensity);
symm_dna_mn_intensity=nanmean(symm_dna_intensity);
symm_mukbef_md_intensity=nanmedian(symm_mukbef_intensity);
symm_mukbef_mn_intensity=nanmean(symm_mukbef_intensity);


%% save data per cell:
ColNames=[{'experiment'}, {'cellname'} , {'symmetry_DNA_by_area'}, {'symmetry_MukB_by_area'} ,{'symmetry_DNA_by_intensity'}, {'symmetry_MukB_by_intensity'},{'rot.axis length'},{'length>5?'},{'area dna>mukb?'}];
xlswrite(strcat(initval.pth_repli,initval.DirSep,initval.expi,'_A062_symmetries.xlsx'),ColNames,'Per Cell','A1');
xlswrite(strcat(initval.pth_repli,initval.DirSep,initval.expi,'_A062_symmetries.xlsx'),expoutnames,'Per Cell','A2');
xlswrite(strcat(initval.pth_repli,initval.DirSep,initval.expi,'_A062_symmetries.xlsx'),celloutnames,'Per Cell','B2');
xlswrite(strcat(initval.pth_repli,initval.DirSep,initval.expi,'_A062_symmetries.xlsx'),symmetries,'Per Cell','C2'); 

%% save summary
RowNames=[{'number of cells'}, ...
          {'DNA_by area, mean'}, {'DNA_by area, median'}, {'MukB_by area, mean'} ,{'MukB_by area, median'} ...
          {'DNA_by intensity, mean'}, {'DNA_by intensity, median'}, {'MukB_by intensity, mean'} ,{'MukB_by intensity, median'}]';
rowout=[
    N_ok;
    symm_dna_mn_area;
    symm_dna_md_area;
    symm_mukbef_mn_area;
    symm_mukbef_md_area;
    symm_dna_mn_intensity;
    symm_dna_md_intensity;
    symm_mukbef_mn_intensity;
    symm_mukbef_md_intensity];
xlswrite(strcat(initval.pth_repli,initval.DirSep,initval.expi,'_A062_symmetries.xlsx'),RowNames,'Summary','A1');
xlswrite(strcat(initval.pth_repli,initval.DirSep,initval.expi,'_A062_symmetries.xlsx'),rowout,'Summary','B1'); 


%figure(4);
%subplot(1,2,1);
%violindata=[symm_dna_area symm_mukbef_area];

%[h,L,MX,MED]=violin(violindata, 'xlabel',{'DNA','MukBef'});
%ylabel('symmetry','FontSize',14)

% subplot(1,2,2);
% plot(symm_dna_area,symm_mukbef_area,'ro'); hold on;
% plot([-1 1], [0 0],'k-');
% plot([0 0],[-1 1] ,'k-');
% xlabel('symmetry ratio DNA, a.u.');
% ylabel('symmetry ratio MukB, a.u.');
% axis equal
% xlim([0 1]);
% ylim([0 1]);
% saveas(gcf,strcat(initval.pth_repli,initval.DirSep,initval.expi,'_A062_symmetries.jpg'));

           