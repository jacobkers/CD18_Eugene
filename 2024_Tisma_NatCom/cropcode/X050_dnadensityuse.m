function X050_dnadensityuse(pathinfo, dryrun);
%JWJK_A:-------------------------------------------------------------------
%Title: 'X050_dnadensityuse'
%Summary:  This code re-saves data as input for the 'Donut code'
%Input: all info from Grid_cell_001t06xy1.mat etc
%Output: same naming, only pictures; excel list of names (string)
%References: code by F. Wu and X. Zheng, reedited by JK'17
%:JWJK_A-------------------------------------------------------------------
if nargin<1, pathinfo = X000_setpath4snapshots,end

if 0
    if isdir(pathinfo.dirdensityanalysis) 
        disp('deleting.....');
        rmdir(pathinfo.dirdensityanalysis,'s');  end
        disp('done deleting');
    
end
if ~ isdir(pathinfo.dirdensityanalysis), 
    mkdir(pathinfo.dirdensityanalysis); 
end
cd(pathinfo.dircrop);
load('cell_namelist.mat','cellnames');
N_cells=length(cellnames);

cd(pathinfo.dirgrid); % folder with selected cells
file = dir('Grid*');

for ci = 1:N_cells;
    cd(pathinfo.dircrop); % where the cropped cells are saved in mat files
    thiscellname=char(cellnames{ci});
    disp(['X050_',num2str(ci),'_working:', thiscellname,'_with', num2str(N_cells-ci), 'to go'])
    ni=strfind(thiscellname, 'chan');
    namepartA=thiscellname(1:ni-1);
    namepartB=thiscellname(ni+4:end); 
    if ~dryrun
        load(['ma' namepartB '.mat']); %cell boundary filenaming
        load(['c1' namepartB '.mat']); %c1 phase image name
        load(['c2' namepartB '.mat']);  %c2 hu 
        load(['c3' namepartB '.mat']);  %c3 ori
        load(['c4' namepartB '.mat']);  %c4 ter
        load(['c5' namepartB '.mat']);  %c5 mukbef optional
        load(['rc2' namepartB '.mat']);  %c2 hu 
        load(['rc3' namepartB '.mat']);  %c3 ori
        load(['rc4' namepartB '.mat']);  %c4 ter
        load(['rc5' namepartB '.mat']);  %c5 mukbef optional
        load(['xc2' namepartB '.mat']);  %c2 hu 
        load(['xc3' namepartB '.mat']);  %c3 ori
        load(['xc4' namepartB '.mat']);  %c4 ter
        load(['xc5' namepartB '.mat']);  %c5 mukbef optional

        if 0
            subplot(2,3,1); pcolor(cellmask); shading flat, colormap bone; axis equal, axis tight, axis off; 
            subplot(2,3,2); pcolor(cellc1); shading flat, colormap bone; axis equal, axis tight, axis off;                        
            subplot(2,3,3); pcolor(cellc2(:,:,2)); shading flat, colormap bone; axis equal, axis tight, axis off;
            subplot(2,3,4); pcolor(cellc3(:,:,2)); shading flat, colormap bone; axis equal, axis tight, axis off;
            subplot(2,3,5); pcolor(cellc4(:,:,2)); shading flat, colormap bone; axis equal, axis tight, axis off;
            subplot(2,3,6); pcolor(cellc5(:,:,2)); shading flat, colormap bone; axis equal, axis tight, axis off;  
            [~]=ginput(1);
        end


        cd(pathinfo.dirdensityanalysis); % folder where to save the data for donut code
        save(['ma' namepartB '.mat'],'cellmask');
        save(['c1' namepartB '.mat'],'cellc1');
        save(['c2' namepartB '.mat'],'cellc2','cellc2_rw', 'cellc2_xt');
        save(['c3' namepartB '.mat'],'cellc3','cellc3_rw', 'cellc3_xt');
        save(['c4' namepartB '.mat'],'cellc4','cellc4_rw', 'cellc4_xt');
        save(['c5' namepartB '.mat'],'cellc5','cellc5_rw', 'cellc5_xt');


        savefullstack=0;
        if savefullstack        
            save(['c2_' cellname '.mat'],'cellc2_fs_raw', 'cellc2_fs_decon','-append');  %optional
        end
    end
end
cd(pathinfo.dircode);    