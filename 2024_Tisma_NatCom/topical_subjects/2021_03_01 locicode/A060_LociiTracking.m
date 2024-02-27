function A060_LociiTracking(initval)
%JWJK_A:-------------------------------------------------------------------
%Building of  locii traces
%
%Summary: code loads image stack per cropped cell and performs spot
%analysis on each.
%
%Approach: first, stack is loaded of rfp or cfp. average of stack is used
%for defining a fixed, smaller ROI. This is then used for 2D-Xcor tracking
%
%Input: movie stack onf cropped cell images; all channels 
%
%Output: collection of traces per cell
%
%References: written by JWJK
%:JWJK_A-------------------------------------------------------------------
sho=1;
close all;
    curpath=pwd;
    txtoutdir=strcat(initval.imagepath,'A60_LociiTraces\');
    %if isdir(txtoutdir), rmdir(txtoutdir,'s');  end
    mkdir(txtoutdir);
    %build list, based on bf images
    %1) count number of cells based on firstframes
    close all;
    codepth=pwd;
    indir =strcat(initval.imagepath,'A50_Cropped\');
    channellabels=initval.SaveLabels;
    %build list, based on bf images
    %1) count number of cells based on firstframes
    cd(indir);
    cell_list = dir('*bf.tif');
    [CellsNo,~]=size(cell_list);               %number of cells
    FramesNo = length( imfinfo('nr001bf.tif'));     %frames
    cd(codepth);
    set(figure(3), 'visible', 'off');
    for ii=1:CellsNo , %ii %for all cells; 
        cell_index=cell_list(ii).name(3:end-6);
        disp(strcat(initval.imagepath,'_Locii Tracking','Cell',cell_index));
        %load movie
        %label
        namebase=[cell_list(ii).name(1:5)  channellabels{5}];
        filname=[namebase '.tif'];  %chromosome
        label_im=double(imread(strcat(indir,filname)));  

        %remove other label areas
        [rr,cc]=size(label_im); xc=round(cc/2); yc=round(rr/2);
        labelval=label_im(yc,xc);    
        label_im(label_im~=labelval)=0; 
        pic0=1.0*(label_im-min(label_im(:)))/(max(label_im(:))-min(label_im(:)));
        Msk0=1.0*bwmorph(pic0, 'dilate',3);
        Msk1=matrix_blur_jk(Msk0,3); %soft edge blur
        [rr,cc]=size(pic0);
         raw_xtraces=[];
         raw_ytraces=[]; 
        for ch=2:length(channellabels)-1 %for all color channels (2,3,4)
            %% obtain the stack
            this_color=channellabels{ch};
                namebase=[cell_list(ii).name(1:5)  this_color];
                filname=[namebase '.tif'];  %chromosome
                for jj=1:FramesNo 
                    CellFrameIm=double(imread(strcat(indir,filname),'Index',jj));
                    if jj==1             
                        ImStack=zeros(rr,cc,FramesNo);               
                    end
                    ImStack(1:rr,1:cc,jj)=Msk1.*double(CellFrameIm);
                end
            %%---------------------
            AvSpotIm=mean(ImStack,3);
            
            [xm,ym]=TrackXYByMax(AvSpotIm); % 1) track center by max to cut fixed sub-roi
            if 0 %sho
                subplot(2,1,1);
                pcolor(AvSpotIm); shading flat; axis equal; axis tight; hold on;
                title('average image: pre-locate spots');
            end
            xtrace=zeros(FramesNo,1);
            ytrace=zeros(FramesNo,1); 
            reftrace_x=zeros(FramesNo,1);
            reftrace_y=zeros(FramesNo,1);
            ok2track=Check_image(AvSpotIm,xm,ym);
            if 1 %ok2track          
                for jj=1:FramesNo      %for all frames
                    FrameIm=ImStack(:,:,jj);
                    [roi,xoff,yoff]=Cut_roi(FrameIm,xm,ym); 
                    switch this_color
                        case 'color1', 
                            %[x1,y1,~,~]=TrackXY_by_2DXCor(roi);
                            [x1,y1,~,~]=TrackXY_by_COM_2Dmoment(roi);
                        case 'color2', 
                            %[x1,y1,~,~]=TrackXY_by_2DXCor(roi);
                            [xa,ya,~,~]=TrackXY_by_COM_2Dmoment(roi);
                            [x1,y1,~,~,~,~,~]=track_xy_by_com_clipmask(roi,xa,ya);
                        case 'color3', 
                             [xa,ya,~,~]=TrackXY_by_COM_2Dmoment(roi);
                            [x1,y1,~,~,~,~,~]=track_xy_by_com_clipmask(roi,xa,ya);
                    end
                    xtrace(jj)=x1;
                    ytrace(jj)=y1;                   
                end
                if 0 %sho
                    subplot(2,1,2);
                    pcolor(roi); shading flat; axis equal; axis tight; hold on;
                    plot(xtrace+0.5,ytrace+0.5,'k-'); hold on;
                    text(1,-2,filname);
                    title('small ROI for tracking')
                    pause(0.2);
                    [~]=ginput(1);
                end
                if ch==2 
                    reftrace_x=xtrace-xtrace(1);
                    reftrace_y=ytrace-ytrace(1); 
                end
                
                raw_xtraces=[raw_xtraces xtrace-xtrace(1)];
                raw_ytraces=[raw_ytraces ytrace-ytrace(1)];
                
                %remove offset;  
                xtrace=xtrace-xtrace(1);
                ytrace=ytrace-ytrace(1);
                 
                %remove drift; save for MSD analysis 
                xtrace_nrm=xtrace-reftrace_x;
                ytrace_nrm=ytrace-reftrace_y;                    
                TraceName=strcat('TraceOfCellNr',cell_index,...
                                      '_',char(channellabels(ch)));                                
                dlmwrite(strcat(txtoutdir,TraceName,'.txt'),[xtrace_nrm ytrace_nrm]);
                if ch==length(channellabels)-1;
                        save(strcat(txtoutdir,'TracesOfCellNr',cell_index,'.mat'),'raw_xtraces','raw_ytraces');
                end
            end
            figure(3);
            subplot(2,1,1);        plot(xtrace); hold on;
            ylabel('X-position');
            title(['A060-Locii Tracking','Cell',num2str(ii)]);
            subplot(2,1,2);        plot(ytrace); hold on;
            ylabel('Y-position');
            xlabel('Time, frames');              
        end 
        pause(0.1);
        legend(channellabels(2:4));
        PlotName=strcat('TraceOfCellNr',num2str(ii,'%03i'), '.jpg');
        saveas(gcf, strcat(txtoutdir,PlotName));        
        %[~]=ginput(1);
        close(gcf);
    end
    cd(curpath)

function okall=Check_image(Im,xm,ym);
% JWKC:
%Check if spot image is regular, if not, give error flag
% :JWKC
Im=Im-median(Im(:));      %subtract background
Im=Im/sum(Im(:))*100;  %normalize
[rr,cc]=size(Im);

edg=6; %spotnot too close to edges?
okedge=(xm>edg)&(xm<(cc-edg))&(ym>edg)&(ym<(rr-edg));
okspot=1;
if okedge %spot contains most of intensity?
    psf2=5;
    spotroi=Im(ym-psf2:ym+psf2,xm-psf2:xm+psf2);
    spotsum=sum(spotroi(:));
    okspot=spotsum>25;
end    
okall=okedge&okspot;


