function verdict=user_judge_cell(extra_im,plot_chx,plot_chy,plot_cellx,plot_celly,counter);
        figure(2);
        if ~isempty(extra_im);
            subplot(1,2,2); pcolor(extra_im); shading flat, axis equal, colormap hot; axis off
            subplot(1,2,1);
        end
        corners=[min(plot_cellx) min(plot_celly);
        max(plot_cellx) min(plot_celly);
        max(plot_cellx) max(plot_celly);
        min(plot_cellx) max(plot_celly)];
        plot(plot_chx,plot_chy,'r-', 'Linewidth',5); hold on;
        plot(plot_cellx,plot_celly,'k--'); hold on;
        xlabel('contour circularity, a.u.');
        ylabel('contour symmetry, a.u.');
        hold on
        axis equal
        axis off
        title(['cell:',num2str(counter),'-what is it?']);

        clr='k';
        shft=3;
        text(corners(1,1),corners(1,2), '1-crescent', 'Color', clr);
        text(corners(2,1),corners(2,2),'2-compact', 'Color', clr);
        text(corners(3,1),corners(3,2),'3-pancake', 'Color', clr);
        text(corners(4,1),corners(4,2),'4-multilobe', 'Color', clr); 
        pause(0.01);
        tic
        [xn,yn,but]=ginput(1);
        clicktime=toc;
        dist=((corners(:,1)-xn).^2+(corners(:,2)-yn).^2).^0.5;
        [mindist,verdict]=min(dist);              
        close(gcf);
        %append to result
