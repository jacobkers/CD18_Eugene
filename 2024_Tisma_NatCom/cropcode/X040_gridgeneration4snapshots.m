function X040_gridgeneration4snapshots(pathinfo)
%JWJK_A:-------------------------------------------------------------------
%Title: X040_gridgeneration4snapshots
%Summary: create summary pictures of cell data
%Input: files 'c2_cell_im00014t04xy1.mat' per cell area, time frame and imaging position; 
%Output: 'Grid_cell_001t01xy1' color tiffs
%References: code by F. Wu and X. Zheng, reedited by JK'17
%:JWJK_A-------------------------------------------------------------------

if nargin<1, pathinfo = X000_setpath4snapshots,end;

%if isdir(pathinfo.dirgrid), rmdir(pathinfo.dirgrid,'s');  end
if isdir(pathinfo.dirmeasurement), rmdir(pathinfo.dirmeasurement,'s');  end

mkdir(pathinfo.dirgrid); 
mkdir(pathinfo.dirmeasurement); 

% define a few parameters
% pxval=15;
% numfr=3;
% midfr=2;

% go to cropped cell folder

cd(pathinfo.dircrop);
load('cell_namelist.mat','cellnames');
% h = waitbar(0,'Please wait...');
N_cells=length(cellnames);
cellindex = [];chrindex = []; oriindex = [];terindex=[];cellsource=[];cellchr=[];chrindexterincluded=[];
for ci = 1:N_cells
        %waitbar(f/steps);
        thiscellname=char(cellnames{ci});

        ni=strfind(thiscellname, 'chan');
        namepartA=thiscellname(1:ni-1);
        namepartB=thiscellname(ni+4:end);        
               
        disp([num2str(ci),'_working:', thiscellname,'_with', num2str(N_cells-ci), 'to go']);
        cd(pathinfo.dircrop);

            % phase contrast images
            load(['c1' namepartB '.mat'],'cellc1');
            phcrop = cellc1;
            B1=phcrop-min(phcrop(:));
            B1=(B1./(max(phcrop(:))-min(phcrop(:))))*255;
            B1R=uint8(cat(3,B1,B1,B1));
            
            % cell mask + measurement
            load(['ma' namepartB '.mat'],'cellmask');
            labcrop = cellmask;
            B2=(labcrop==1);
            Ilab=label(logical(labcrop),2,100,10000);
            if max(Ilab(:))>0
                Ilab1 = (Ilab==1);
                % measure the mask
                data=measure(Ilab1,[],{'size','feret','center','minimum','maximum'});
                celldat=cat(2,str2double(thiscellname),double(data.size)',double(data.feret)',(data.center)',(data.minimum)',(data.maximum)');
                %overlay mask+phase contrast
%                 if 1 %fi==61
%                     dipshow(Ilab1); hold on;
%                     plot(data.Center(1),data.Center(2)); hold off;
%                     fi
% %                     size(Ilab1)
% %                     data.Center
%                     pause(0.7);
%                  end
                B2b=B2.*255;
                B2R=uint8(cat(3,B2b,B2b,B2b));
                [croph, cropw]=size(labcrop);
                B3R=uint8(cat(3, B2b, zeros(croph,cropw), B1)); % RGB: blue phase image with red identified cell
                
                % Take the chromosome image and measure the size
                load(['c2' namepartB '.mat'],'cellc2'); 
                hucrop = cellc2;
                B4=hucrop(:,:,2);  %assuming three planes[JK]               
                B4b=B4-min(B4(:));
                B4c=(B4b./(max(B4(:))-min(B4(:))))*255;
                B4R=uint8(cat(3,B4c,B4c,B4c));
                B5=[];
                [~,B5c,~]=find_obj2(B4,8,B2,1,3,0.1);
                B5=double(B5c);
                [chrdat1,Chrbi,Chrsum,Chrcent] = find_chr1(B5); 
                if ~isempty(chrdat1)
                    numchr=size(chrdat1.size,2);
                    chrdat=[double(ones(numchr,1))*str2double([thiscellname]) double(Chrsum) double(chrdat1.feret)' double(Chrcent)];
                else
                    numchr=1;
                    chrdat=[str2double([thiscellname]) 0 0 0 0 0 0 0 0];
                end
                B5d=double(Chrbi)*255;
                B5R=uint8(cat(3,B5d,B5d,B5d));            
                B6b=double(B5d);
                B6R=uint8(cat(3,B6b*255,zeros(croph,cropw),B2b));
                
                % Take the Ori/Ter data
                load(['c3' namepartB '.mat'],'cellc3');
                load(['c4' namepartB '.mat'],'cellc4');              
                oricrop = cellc3;
                tercrop = cellc4;
                B70=max(oricrop,[],3);
                B7=B70;
                B7b=(B7-min(B7(:))).*255/(max(B7(:))-min(B7(:)));
                [B7d,~,oridat1]=find_obj2(B7,8,B2,1,3,0.2);
                NumCellLabel=BuildNumericCellLabel(thiscellname);
                
                %numori=size(oridat1.center,2);
                if ~isempty(oridat1)
                    numori=size(oridat1.center,2);
                    oridat=[double(ones(numori,1))*NumCellLabel double(oridat1.center)'];
                else
                    numori=0;
                    oridat=[NumCellLabel 0 0];
                end
                %oridat=[double(ones(numori,1)) double(oridat1.center)'];
                B7R=uint8(cat(3,B7b,B7b,B7b));
                B80=max(tercrop,[],3);
                B8=B80(:,:,1);
                B8b=(B8-min(B8(:))).*255/(max(B8(:))-min(B8(:)));
                
                [B8d,B8c,terdat1]=find_obj2(B8,8,B2,1,3,0.2);
                if ~isempty(terdat1)
                    numter=size(terdat1.center,2);
                    terdat=[double(ones(numter,1))*NumCellLabel double(terdat1.center)'];
                else
                    numter=0;
                    
                    terdat=[NumCellLabel 0 0];
                end
               
                B8R=uint8(cat(3,B8b,B8b,B8b));
                B9R=cat(3,uint8(B7d*255),uint8(B8d*255),uint8(B2b));
                
                % Only save the results for single chromosomal cells
                if numori > 0 && numter > 0
                    % What if we put terlabel into hulabel 
                    % terlabel B8d hulabel Chrbi
                    obchromo = B8d + Chrbi;
                    obchromo = obchromo >= 1;
                    [Chrdatnew,Chrbinew,Chrsumnew,Chrcentnew] = find_chr1(obchromo);
                    chrdatterincluded=[double(ones(numchr,1))*str2double(thiscellname) double(Chrsumnew) double(Chrdatnew.feret)' double(Chrcentnew)];
                    B5dnew=double(Chrbinew)*255;
                    B12R=uint8(cat(3,B5dnew,B5dnew,B5dnew)); % chromosome + ter binary
                    B11R=uint8(cat(3,uint8(B8d*255),uint8(Chrbi*255),B2b)); % overlay chromosome and ter
                    B10R = uint8(cat(3,uint8(B7d*255),uint8(Chrbi*255),B2b)); % overlay chromosome and ori
                                       
                    % The grid
                    Grid=uint8(cat(1,cat(2,B1R,B2R,B3R),cat(2,B4R,B5R,B6R),cat(2,B7R,B8R,B9R),cat(2,B10R,B11R,B12R)));
                    
                    cd(pathinfo.dirgrid);
                    imwrite(Grid,cat(2,'Grid_',thiscellname,'.tif'),'tiff','Compression','none');

                    % The measurement results to be saved
                    cellindex = cat(1,cellindex,celldat);
                    chrindex = cat(1,chrindex,chrdat);
                    terindex=cat(1,terindex,terdat);
                    oriindex=cat(1,oriindex,oridat);
                    chrindexterincluded=cat(1,chrindexterincluded,chrdatterincluded);
                    cellchr=cat(1,cellchr,[str2double(thiscellname) celldat chrdat(2:7) chrindexterincluded(2:7)]);                     
                end  
            end
end
if ~isdir(pathinfo.dirmeasurement), mkdir(pathinfo.dirmeasurement); end
cd(pathinfo.dirmeasurement);
save(cat(2,['alldata_' pathinfo.date '.mat']),'cellindex','chrindex','chrindexterincluded','oriindex','terindex','cellchr');   
cd(pathinfo.dircode);             
          
            
        
