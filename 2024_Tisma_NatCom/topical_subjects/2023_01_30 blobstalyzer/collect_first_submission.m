 function collect_first_submission(users)
 runID=['_run_', datestr(now,'yyyy_mm_dd_HH_MM')];
    %% 1 run per image
    for uzi=1:length(users)
        figure('Units','normalized','Position',[0 0 0.5 0.9]);
        usr=char(users{uzi});
        cd(usr);
        click_list=dir('blobstalysis*.mat');
        allblobclass=[];
        for ii =1:length(click_list);    
                load(click_list(ii).name,'typecodes','blobclass','headers');
                allblobclass=[allblobclass ;...
                                blobclass];
        end
        image_ID=allblobclass(:,3);
        blob_ID=allblobclass(:,2);
        binax=[1:1:4];
        BlobLegend=([{'1.donut'},{'2.crescent'},{'3.compact'},{'4.other'}]);
        for jj=1:4
            sel=find(image_ID==jj);
            thisimage_blob_ID=blob_ID(sel);
            thisimage_hist=hist(thisimage_blob_ID,binax);
            N_total=sum(thisimage_hist);
            thisimage_hist=thisimage_hist/N_total*100;
            subplot(4,1,jj);
            for ki=1:4
                bar(binax(ki),thisimage_hist(ki)); hold on;
                title([usr, ':image:', num2str(jj), '; N=',num2str(N_total)]);
            end
            xlabel('type');
            ylabel('occurence, %');
            legend(BlobLegend, 'Location', 'EastOutside');
        end
        cd ..;
        saveas(gcf,[usr '_clicks.jpg']);
        plot2svg([usr '_click.svg'], gcf);
        ColNames=[{'index'} , {'type'}, {'source image'} ,{'blob index'}, {'clicktime'}];
        xlswrite([usr '_clicks_1st submission',runID,'.xlsx'],ColNames,'Sheet1','A1');
        xlswrite([usr '_clicks_1st submission',runID,'.xlsx'],allblobclass,'Sheet1','A2');

        if 0
        % just for fun
        figure(length(users)+1);
        subplot(1,2,uzi);
        binax=0:0.15:5;
        clicktimes=allblobclass(:,5);
        clicktimehist=hist(clicktimes,binax); hold on
        bar(binax, clicktimehist);
        legend(usr);
        title('clicktimes');
        %subplot(2,1,2);
        %plot(clicktimes, 'o-'); hold on;
        %legend(users);
        saveas(gcf,[usr '_click_times.jpg']);
        end
    end