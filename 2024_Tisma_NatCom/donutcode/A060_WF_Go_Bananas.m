function A060_WF_Go_Bananas(batchrunindex)
%JWJK_B:-------------------------------------------------------------------
% checking crescent shape of cells. 
%following 'D:\jkerssemakers\Dropbox\CD_recent\BN_CD22_Tisma\analysis\CD20_Cells\topical_projects\2022_11_29 cell contours'
%
%:JWJK_B-------------------------------------------------------------------

if nargin<1, batchrunindex=-101.2; end
initval=A000__WF_Get_JacobPathsandExperiments(batchrunindex);

%actions.re_load=0; %for re-loading oufti-data
actions.re_run=0;  %for re-doing analysis (otherwise, just post-plotting existing)


%for user-click classification
actions.user_judge=0;       %
actions.re_judge_append=0;  %only when user_judge
%auto-overview:
actions.plot_shapes=1;

%set and transfer paths
curpth=pwd; 
cd ..; addpath(genpath('topical_projects\2022_11_29 cell contours\')); cd(curpth);

%if isdir(imoutdir), rmdir(imoutdir,'s');  end
banana_init.exp_label=initval.expi;
banana_dir=strcat(initval.resultpath,'A060_CellImages_Bananalyzer\',initval.DirSep);
if ~isdir(banana_dir), mkdir(banana_dir);  end
banana_init.datapath_in=banana_dir;
banana_init.datapath_out=banana_dir;

%settings
banana_init.blowup=10;  %over-sample images for smoother contours
banana_init.total_sampling =1400;  %to do all, use number> number of cells
banana_init.user_sampling=1400;
banana_init.px2um = initval.nmperpixel/1000; 
    %% Microscope resolution (um to pixels conversion) - modify this value when required
    %px2um = 0.07;  Dumbledore 1.5X pixel size
    %px2um = 0.11;  Dumbledore pixel size
    %px2um = 0.04;  Snape 1.5X pixel size
    %px2um = 0.065; Snape pixel size

plotinfo.non_overlap=0.015; %in units of symmetry params
plotinfo.scalar=2000; %the smaller, the larger the cells in the plot...  
    
    
close all; pause(0.1);
if actions.re_run   
initval.Cell_Labels;
allframes=length(initval.Cell_Labels);
Nframes=min([200E6 allframes]);
goodcount=0;
shape_data=[];

plotinfo.used_xy=[];
counters.gen = 0;
if actions.user_judge
%check existing data:
    if  actions.re_judge_append  %add clicks to existing data   
        load([banana_init.datapath_out, banana_init.exp_label, '_shape_data.mat']);
        clicked=find(shape_data(:,3)>0);
        counters.last_click_count=length(shape_data(clicked,1));
        counters.click=counters.last_click_count;
    else %start clicking all over
        shape_data=[];
        counters.click=0;
        counters.last_click_count=0;  
    end  
else %just re-run, no new clicks or saves
    load([banana_init.datapath_out, banana_init.exp_label, '_shape_data.mat']);
end


for jj=1:Nframes
    cellno=char(initval.Cell_Labels{jj});
     disp(strcat('assessing bananaism:cell',num2str(jj),'and',num2str(Nframes-jj+1),'cells to go:'));
    %get data
     CellName=strcat('ResultsOfCell',cellno,'.mat'); 
     MatFilePath=strcat(initval.resultpath,'ResultsPerCellMatlab',initval.DirSep);
      load(strcat(MatFilePath,CellName));               
      if (GeneralCellProps.Okayproduct|initval.PassAllCells)&Aligned.BP.CorrectionOK
          
          counters.gen = counters.gen + 1;
          %BW_cell=cell_mask;
          work_im_cell_BW=blow_up_image_or_xy(cell_mask, banana_init.blowup, 1);
          
          [cellX,cellY]=BW_to_contour(work_im_cell_BW);
          %prepare atresholded chroosome picture:
          
          
          work_im_Chro=blow_up_image_or_xy(Chromosome.picture, banana_init.blowup, 0);
          BW_chro=get_single_shape(work_im_Chro);
          
          
          %get and clean edges (see
          %'clusters')--------------------------------
        
        bwstruct=bwconncomp(BW_chro,8);    %finds 8-fold connected regions.
        area_props=regionprops(bwstruct,...
        'Centroid', 'Area','Perimeter', 'Orientation', 'Circularity');
  
        %HANDLE MULTIPLE OBJECTS.....possible average x0?
        N_areas=length(area_props);
        for ai=1:1 %N_areas
            %get clean_sorted contours: 
            [chrX,chrY]=BW_to_contour(BW_chro,'nearest',1);
            %prepare:
            extra_im=work_im_Chro;
            [shape_data,counters,plotinfo]=get_symmetry_and_user_info(chrX,chrY,cellX,cellY,BW_chro,extra_im, banana_init,actions,shape_data,counters,plotinfo);           
        end
         pause(0.01);        
         dum=1;
          
      end
      if actions.user_judge  %save only when re-assessing clicks
          %note that data is re-saved after every frame!
            save([banana_init.datapath_out, banana_init.exp_label, '_shape_data.mat'], 'shape_data');
      end
      
end
end
%post-processing section:
AT030_Bananalyzer_Post_Process(banana_init);


function BW_chro=get_single_shape(work_im_Chro);
%treshold image, if multiple shapes emerge, lower treshold

[~,thr_min]=Find_treshold_MD_V2020(work_im_Chro,0);
mx=max(work_im_Chro(:));
new_scalar=0.5;

multiplets=1;
while multiplets
    current_thr=new_scalar*(mx-thr_min)+thr_min;
    BW_chro=work_im_Chro>current_thr;
    bwstruct=bwconncomp(BW_chro,8);    %finds 8-fold connected regions.
    area_props=regionprops(bwstruct,...
    'Centroid');
%decrease treshold
    if length(area_props)>1 & current_thr>thr_min  %reject multiple objects
            new_scalar=new_scalar-0.1;
    else
            multiplets=0;
    end
 end


