 function A050_Cellcropper_Loci(initval);
%JWJK_A:-------------------------------------------------------------------
%Split microscope data in movies cropped per cell
%
%Summary: Crop pre-determined cell areas and save as separately labeled images
%
%Approach: program loads pre-determined cell data; cropps cell areas
%accordingly and saves separate per color
%
%Input: pre-detrmined cell screening data
%
%Output: cell images
%
%References: written by F.Wu/X. Zheng. Edited &annotated by JWJK
%
%:JWJK_A-------------------------------------------------------------------
%dip init
if 1
addpath('C:\Program Files\DIPimage 2.9\common\dipimage');
dip_initialise;
dipsetpref('ImageFilePath', 'C:\Program Files\DIPimage 2.8\images');
end



codepth=pwd;
    imindir =initval.imagepath;  %here it finds the files
    imoutdir=strcat(imindir,'A50_Cropped\');
     if isdir(imoutdir), rmdir(imoutdir,'s');  end
     mkdir(imoutdir);
    load([initval.matfilesinoutdir, 'A040_BN_labels.mat']);
    
    color1_exists=(length(frsloci_color1)>0);
    color2_exists=(length(frsloci_color2)>0);
    color3_exists=(length(frsloci_color3)>0);
    
    CS=initval.cropping_edge-1; %cropping edge
    Nframes=min([length(frs), initval.shortset]);   
    for f=1:Nframes;
        disp(strcat(initval.expname,'_',initval.subdir,'_Cropping From Movieframe',num2str(f)));
        cd(imindir);
        I=imread(frs(f).name);              %phase images
        if color1_exists, J=imread(frsloci_color1(f).name);  end    
        if color2_exists, L=imread(frsloci_color2(f).name); end      
        if color3_exists, K=imread(frsloci_color3(f).name); end
        %do the cropping, CS=20 pts around area
        cd(imoutdir);                
        for c=1:size(celldat1,1);
            namebase=cat(2,'Frame',num2str(f),'_','nr',num2str(celldat1(c,1),'%03i'));  
            namebase=cat(2,'nr',num2str(celldat1(c,1),'%03i'));
            %set space around cell area of interest
            [rr,cc]=size(I);
            loy=max([1  celldat1(c,11)-CS]);
            hiy=min([rr celldat1(c,13)+CS]);
            lox=max([1  celldat1(c,10)-CS]);
            hix=min([cc celldat1(c,12)+CS]);
                                           
            if f==1;  %label image (averaged mask)              
                Hcrop=uint16(lab1(loy:hiy,lox:hix));
                imwrite(Hcrop,cat(2,namebase,'lab','.tif'),'tiff','Compression','none','WriteMode','Append');
            end                  
            %brightfield or phase image
            Icrop=I(loy:hiy,lox:hix);
            imwrite(Icrop,cat(2,namebase,'bf','.tif'),'tiff','Compression','none','WriteMode','Append');
            %color 1
            if color1_exists
                Jcrop=J(loy:hiy,lox:hix);
                imwrite(Jcrop,cat(2,namebase,'color1','.tif'),'tiff','Compression','none','WriteMode','Append');
            end
            %color 2 
            if color2_exists
                Lcrop=L(loy:hiy,lox:hix);
                imwrite(Lcrop,cat(2,namebase,'color2','.tif'),'tiff','Compression','none','WriteMode','Append');
            end
            %color 2 
            if color3_exists
                Kcrop=K(loy:hiy,lox:hix);
                imwrite(Kcrop,cat(2,namebase,'color3','.tif'),'tiff','Compression','none','WriteMode','Append');
            end
                      
        end
    end
    cd(initval.codepth);
disp('done')
