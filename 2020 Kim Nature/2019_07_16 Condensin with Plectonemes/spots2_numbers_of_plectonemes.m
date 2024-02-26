function spots2_numbers_of_plectonemes(info_DNA,info_Cnd,kymo_DNA,kymo_Cnd,SaveName);
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
    plectonemic_content=zeros(FF,1);
    framespan=5;
    for ii=1:FF
        sel=find((info_DNA.pos_frameno>=ii-framespan)&...
            (info_DNA.pos_frameno<=ii+framespan));
        if~isempty(sel);
            plectonemes_per_frame(ii)=length(sel)/(2*framespan+1);
            plectonemic_content(ii)=sum(info_DNA.content_perspot_meas(sel))/(2*framespan+1);
        end
    end
    
    subplot(1,3,1); pcolor(kymo_DNA); shading flat; colormap hot; 
        title('DNA'); 
    subplot(2,3,3); plot(plectonemes_per_frame,'-');
        title('number of loops/plectonemes');
        ylabel('plectoneme #');
        xlabel('frame no.');
    subplot(2,3,6); plot(plectonemic_content,'r-');
        title('content of loops/plectonemes');
        ylabel('DNA percentage in loops');
        xlabel('frame no.');
        
%         saveas(gcf,strcat(SaveName, '_plectonemecounts.jpg'),'jpg');   
%         close(gcf);
    dum=1;