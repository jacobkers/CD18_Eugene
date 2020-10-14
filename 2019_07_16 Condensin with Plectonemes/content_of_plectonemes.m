function A025_results=content_of_plectonemes(info_DNA,kymo_DNA,SaveName,expinfo);
%JWJK_B:----[add ABCorC*----------------------------------------------------
%Title: %count plectonemes. 
%Summary: %count the number of plectonemes per frame; allow for averaging
%over a 'span'
%Output: plot
%References: CD lab, project Eugene Kim, written by Jacob Kers, 2019
%:JWJK_B-----[add ABCorC*---------------------------------------------------
    if 1
    %some detection diagnosis
        content_est=(info_DNA.content_perspot_meas);
        bins=0:2:100;
        contents=hist(content_est,bins);
        bar(bins,contents);
        dum=1;
    end

    framespan=0;
    mincontent=0;

    %1 number of plectonemes per frame
    FF=max(info_DNA.pos_frameno);
    plectonemes_per_frame=zeros(FF,1);
    plectonemic_content_total=zeros(FF,1);
    plectonemic_content_brightest=zeros(FF,1);
    plectonemic_mean=zeros(FF,1);
    plectonemic_content_sorted=zeros(FF,10);

    
    for ii=1:FF
        sel=find((info_DNA.pos_frameno>=ii-framespan)&...
            (info_DNA.pos_frameno<=ii+framespan)&...
            (info_DNA.content_perspot_meas>mincontent));
        if~isempty(sel);
            plectonemes_per_frame(ii)=round(length(sel)/(2*framespan+1));
            plectonemic_content_total(ii)=sum(info_DNA.content_perspot_meas(sel))/(2*framespan+1);
            plectonemic_content_brightest(ii)=max(info_DNA.content_perspot_meas(sel));
            plectonemic_content_mean(ii)=plectonemic_content_total(ii)/plectonemes_per_frame(ii);
            
            %make sorted map
            sub_fr=1;
            plectonemic_sorted_spanbuf=NaN(1+2*framespan,10);
            for jj=ii-framespan:ii+framespan
                subsel=find((info_DNA.pos_frameno==jj)&...
                    (info_DNA.content_perspot_meas>mincontent));
                if ~isempty(subsel);
                    sortedcontent=sort(info_DNA.content_perspot_meas(subsel),'descend');
                    if length(sortedcontent)>10, %crop
                        sortedcontent=sortedcontent(1:10);
                    end
                        plectonemic_sorted_spanbuf(sub_fr,1:length(sortedcontent))=sortedcontent;
                end
                plectonemic_content_sorted(ii,:)=nanmean(plectonemic_sorted_spanbuf);
                sub_fr=sub_fr+1;
            end
        end
    end
    
    %% save average initial 50 and final 50 frames: DNA content
    A025_results.content_loops_first50=nanmean(plectonemic_content_total(1:50));
    A025_results.content_loops_last50=nanmean(plectonemic_content_total(end-50+1:end));
    A025_results.tetherlength=info_DNA.general_tetherlength;
    
    %% plotting
    [frs,Lx]=size(kymo_DNA);
    subplot(3,2,1); pcolor(kymo_DNA'); shading flat; colormap hot; hold on;
    subplot(3,2,2); 
    pcolor(kymo_DNA'); shading flat; colormap hot; hold on;
    content_sel=find((info_DNA.content_perspot_meas>mincontent));
    plot( info_DNA.pos_frameno(content_sel),info_DNA.pos_X_subpix(content_sel)+expinfo.channelshift, 'bo','Markersize',2); hold on;  

    subplot(3,2,3); 
    plot(plectonemes_per_frame,'k-', 'LineWidth', 2); hold on;
    plot(0*plectonemes_per_frame+1,'k--', 'LineWidth', 0.5);
        ylabel('plectoneme #');
        xlabel('frame no.');
        xlim([1 frs ]);
        ylim([0 6]);
    subplot(3,2,5); 
        %area(plectonemic_content_sorted,'LineStyle','none');
        plot(plectonemic_content_total,'r-'); hold on;
        plot(plectonemic_content_brightest,'k-'); hold on;
        ylabel('cumulative %');
        xlabel('frame no.');
        xlim([1 frs ]);
        ylim([1 120]);
        %colormap jet;
        saveas(gcf,strcat(SaveName, '_plectonemecounts.jpg'),'jpg'); 
        save(strcat(SaveName, '_plectonemecounts.mat'), ...
            'plectonemes_per_frame',...
            'plectonemic_content_total',...
            'plectonemic_content_brightest',...
            'plectonemic_content_sorted',...
            'A025_results');
        close(gcf);
    dum=1;