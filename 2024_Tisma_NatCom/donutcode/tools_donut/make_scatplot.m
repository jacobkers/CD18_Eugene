 function scatplotdata=make_scatplot(map_x,map_y,map_z,labelx,labely,labelz,splitperc,modus);
 %make scatter plot with two classes, color-split by 'splitperc'
 %use map_z as marker size
    [rb,cb]=size(map_x);
    plot(-10,-10,'ro', 'MarkerFaceColor', 'r'); hold on;
    plot(-10,-10,'bo', 'MarkerFaceColor', 'b'); hold on;
    plot([0 100], [0 100], 'k-');
    axis equal; 
    xlim([0 100]);
    ylim([0 100]);
    scatplotdata=[];
    for ii=1:rb
        for jj=1:cb
            scat_x=map_x(ii,jj);
            scat_y=map_y(ii,jj);
            scat_sz=ceil(map_z(ii,jj)/10);
            scatplotdata=[scatplotdata ; [scat_x scat_y scat_sz]];
            if strcmp(modus, 'color')
                if map_z(ii,jj) <splitperc
                    scatcolor='b';
                else
                    scatcolor='r';
                end
                %plot(scat_x,scat_y, [scatcolor, 'o'],'MarkerSize',scat_sz, 'MarkerFaceColor', scatcolor); hold on;
                plot(scat_x,scat_y, [scatcolor, 'o'],'MarkerSize',scat_sz); hold on;
            else
                plot(scat_x,scat_y, 'o', 'MarkerEdgeColor', [0.5 0.5 0.5],'MarkerSize',scat_sz); hold on;
            end
        end
    end
    
    xlabel([labelx '-percentage, %']);
    ylabel([labely '-percentage, %']);
    xlim([0 100]);
    ylim([0 100]);
    %legend([labelz ':+ >',num2str(splitperc),'%'],[labelz ':- <',num2str(splitperc),'%'], 'Location', 'NorthOutside'); 
    axis square;

           