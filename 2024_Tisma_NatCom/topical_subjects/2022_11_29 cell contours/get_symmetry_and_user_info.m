function [shape_data,counters,plotinfo]=get_symmetry_and_user_info(chrX,chrY,cellX,cellY,BW,extra_im, init,actions,shape_data,counters,plotinfo)
%% Use: chr_X,Y; cellX,Y and a chromosome mask to add to exisitng data
%-----------------------------------

bwstruct=bwconncomp(BW,8);    %finds 8-fold connected regions.
area_props=regionprops(bwstruct,...
'Centroid', 'Area','Perimeter', 'Orientation', 'Circularity');
if length(area_props)==1  %reject multiple objects
    x0=area_props.Centroid(1);
    y0=area_props.Centroid(2);
    alpha=area_props.Orientation;
    circularity=area_props.Circularity;           
    
    %rotate long axis horizontal:
    [chrXr,chrYr]=Rotate_Points(x0,y0,chrX,chrY,alpha);
    [cellXr,cellYr]=Rotate_Points(x0,y0,cellX,cellY,alpha);
    %equalize points:
    if 1
        pts=length(chrX);
        [chrXr0,chrYr0, ~]=B002_EqualizeAlongContour(chrXr,chrYr,pts);
        [cellXr0,cellYr0, ~]=B002_EqualizeAlongContour(cellXr,cellYr,pts);
    else
        chrXr0=chrXr;
        chrYr0=chrYr;
        cellXr0=cellXr;
        cellYr0=cellYr;
    end
    %Get middle y (as opposed to ycom
    lox=min(chrXr0);
    hix=max(chrXr0);
    slots=linspace(lox,hix,4);
    y_M=mean(chrYr0(chrXr0>=slots(2)&chrXr0<slots(3)));

    %get symmetry measures:
    ups=length(find(chrYr0-y_M>0));
    downs=length(find(chrYr0-y_M<=0));
    if ups<downs %flip plot!
        flipsign=-1;
    else
        flipsign=1;
    end

    contour_circ=circularity;   
    contour_symm=1-abs(ups-downs)/(ups+downs);
    
   
    
    plotscalar=plotinfo.scalar*init.blowup;  %sets size of cells-to-plot                
    switch 1
        case 1
            plot_chx=(chrXr0-x0)/plotscalar+contour_circ;
            plot_chy=flipsign*(chrYr0-y0)/plotscalar+contour_symm;
            plot_cellx=(cellXr0-x0)/plotscalar+contour_circ;
            plot_celly=flipsign*(cellYr0-y0)/plotscalar+contour_symm;
        case 2
            plot_chx=(chrXr0-x0)/plotscalar+contour_excessL;
            plot_chy=flipsign*(chrYr0-y0)/plotscalar+contour_symm;
            plot_cellx=(cellXr0-x0)/plotscalar+contour_excessL;
            plot_celly=flipsign*(cellYr0-y0)/plotscalar+contour_symm;
    end            

    if actions.user_judge & (counters.gen> counters.last_click_count) & counters.click<init.user_sampling
        %user click section
        verdict=user_judge_cell(extra_im, plot_chx,plot_chy,plot_cellx,plot_celly,counters.gen);   
        counters.click=counters.click+1;
        this_cell_entry=[contour_circ contour_symm verdict];
        shape_data(counters.click,:)=this_cell_entry;
    else  %automatic
        verdict = 0;
        this_cell_entry=[contour_circ contour_symm verdict];
        if ~actions.re_judge_append
            shape_data(counters.gen, :)=this_cell_entry;
        else %add only when new
            if counters.gen>counters.last_click_count
            shape_data(counters.gen, :)=this_cell_entry;
            end
        end
    end  



    %% collective plotting (non-overlap): 
    occupied=check_plotspot(contour_circ,contour_symm,plotinfo.used_xy, plotinfo.non_overlap);
    if actions.plot_shapes & ~occupied
        figure(1);
        verdict_simbol=shape_data(counters.gen, 3);
        switch verdict_simbol
            case 0, clr='r';  lw=1;
            case 1, clr='y-'; lw=2;
            case 2, clr='g-'; lw=1;
            case 3, clr='b-'; lw=1;
            case 4, clr='m-'; lw=1;
        end

        plot(plot_cellx,plot_celly,'k-', 'LineWidth', lw); hold on;
        plot(plot_chx,plot_chy,clr); hold on;
        xlabel('contour circularity, a.u.');
        ylabel('contour symmetry, a.u.');
        hold on
        axis tight
        axis equal
        %[~]=ginput(1);
        %pause(0.1);
        plotinfo.used_xy=[plotinfo.used_xy; [contour_circ contour_symm]];
    end
else
    %still multiple objects, add to DB for proper percentages
    if actions.user_judge & (counters.gen> counters.last_click_count) & counters.click<init.user_sampling
        %user click section 
        counters.click=counters.click+1;
        this_cell_entry=[1.1 1.1 -1];
        shape_data(counters.click,:)=this_cell_entry;
    else  %automatic
        this_cell_entry=[1.1 1.1 -1];
        if ~actions.re_judge_append
            shape_data(counters.gen, :)=this_cell_entry;
        else %add only when new
            if counters.gen>counters.last_click_count
            shape_data(counters.gen, :)=this_cell_entry;
            end
        end
    end  
end