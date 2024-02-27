function A070_LociiShow(initval)
%JWJK_A:-------------------------------------------------------------------
%
%Building of  locii traces
%
%Summary: code loads image stack per cropped cell and performs spot
%analysis on each.
%
%Approach: 
%
%Input: movie stack onf cropped cell images; all channels 
%
%Output: collection of traces per cell
%
%References: written by JWJK
%:JWJK_A-------------------------------------------------------------------
action.savefirstandlastframe=1;
action.showallframes=0;

close all;
set(figure(1), 'visible','off');
codepth=pwd;

trace_indir=strcat(initval.imagepath,'A60_LociiTraces\');
indir =strcat(initval.imagepath,'A50_Cropped\');
picoutdir =strcat(initval.imagepath,'A70_FirstLastFramePics\');
%if isdir(picoutdir ), rmdir(txtoutdir,'s');  end
mkdir(picoutdir );


chans=initval.SaveLabels;
%build list, based on bf images
%1) count number of cells based on firstframes
cd(indir);
cell_list = dir('*bf.tif');
[CellsNo,~]=size(cell_list);               %number of cells
FramesNo = length( imfinfo('nr001bf.tif'));     %frames

cd(codepth);

for ii=1:CellsNo
    %load movie
    %label
    namebase=[cell_list(ii).name(1:5)  chans{5}];
    cell_index=cell_list(ii).name(3:end-6);
    disp(strcat(initval.imagepath,'A70_Locii Showing','Cell',cell_index));
    filname=[namebase '.tif'];  %chromosome
    label_im=double(imread(strcat(indir,filname)));  
   
    tracedataname=['TracesOfCellNr' namebase(3:end-3) '.mat'];
    load([trace_indir,tracedataname], 'raw_xtraces','raw_ytraces');
    
    %remove other label areas
    [rr,cc]=size(label_im); xc=round(cc/2); yc=round(rr/2);
    labelval=label_im(yc,xc);    
    label_im(label_im~=labelval)=0; 
    pic0=1.0*(label_im-min(label_im(:)))/(max(label_im(:))-min(label_im(:)));
    Msk0=1.0*bwmorph(pic0, 'dilate',3);
    Msk0(Msk0==0)=NaN;
    firstlastframe=[];
    for jj=1:FramesNo 
        if action.showallframes || jj==1 |jj==FramesNo
            pic=pic0;
                %colors, masked with label
                for ch=1:4
                    namebase=[cell_list(ii).name(1:5)  chans{ch}];
                    filname=[namebase '.tif'];  %chromosome
                    CellFrameIm=double(imread(strcat(indir,filname), 'Index', jj));  
                    CellFrameIm=Msk0.*(CellFrameIm-min(CellFrameIm(:)));
                    
                    
                    normim=(CellFrameIm-min(CellFrameIm(:)))/(max(CellFrameIm(:))-min(CellFrameIm(:)));
                    pic=[pic normim];
                end
            pcolor(pic); shading flat; colormap bone; axis equal, axis tight
            title([filname num2str(jj)]);
            pause(0.03);
        end
        if (action.savefirstandlastframe &(jj==1|jj==FramesNo))
            firstlastframe=[firstlastframe ; pic];
            if jj==FramesNo                
                picname=[filname(1:end-4) '.jpg'];               
                subplot(2,1,1);
                pcolor( firstlastframe); shading flat; colormap bone; axis equal, axis tight
                title([filname(1:end-4) 'first-last-frame']);
                subplot(2,2,3); plot(initval.pix2nm*raw_xtraces);
                ylabel('x-position,nm');   xlabel('time, frames'); 
                ylim([-60 60]);
                legend([{'color1'},{'color2'},{'color3'}], 'Location', 'SouthOutside');
                subplot(2,2,4); plot(initval.pix2nm*raw_ytraces);
                ylabel('y-position,nm');   xlabel('time, frames'); 
                ylim([-60 60]);
                legend([{'color1'},{'color2'},{'color3'}], 'Location', 'SouthOutside');
                saveas(gcf,[picoutdir,picname]);
            end
        end
    end
end

