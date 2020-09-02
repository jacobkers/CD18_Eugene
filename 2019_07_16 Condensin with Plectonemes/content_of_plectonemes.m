function content_of_plectonemes(info_DNA,kymo_DNA,SaveName);
%JWJK_B:----[add ABCorC*----------------------------------------------------
%Title: %count plectonemes. 
%Summary: %count the number of plectonemes per frame; allow for averaging
%over a 'span'
%Output: plot
%References: CD lab, project Eugene Kim, written by Jacob Kers, 2019
%:JWJK_B-----[add ABCorC*---------------------------------------------------
    
    %1 number of plectonemes per frame
    FF=max(info_DNA.pos_frameno);
    plectonemes_per_frame=zeros(FF,1);
    plectonemic_content_total=zeros(FF,1);
    plectonemic_content_brightest=zeros(FF,1);
    framespan=4;
    for ii=1:FF
        sel=find((info_DNA.pos_frameno>=ii-framespan)&...
            (info_DNA.pos_frameno<=ii+framespan));
        if~isempty(sel);
            plectonemes_per_frame(ii)=round(length(sel)/(2*framespan+1));
            plectonemic_content_total(ii)=sum(info_DNA.content_perspot_meas(sel))/(2*framespan+1);
            plectonemic_content_brightest(ii)=max(info_DNA.content_perspot_meas(sel));
        end
    end
    [frs,Lx]=size(kymo_DNA);
    subplot(2,1,1); pcolor(kymo_DNA'); shading flat; colormap hot; 
        title('DNA'); 
    subplot(2,2,3); plot(plectonemes_per_frame,'-');
        title('number of loops/plectonemes');
        ylabel('plectoneme #');
        xlabel('frame no.');
        xlim([1 frs ]);
        ylim([0 8]);
    subplot(2,2,4); 
        plot(plectonemic_content_total,'r-'); hold on;
        plot(plectonemic_content_brightest,'k-'); hold on;
        title('content of loops/plectonemes');
        ylabel('DNA percentage in loops');
        xlabel('frame no.');
        xlim([1 frs ]);
        ylim([1 100]);
        saveas(gcf,strcat(SaveName, '_plectonemecounts.jpg'),'jpg');   
        close(gcf);
    dum=1;