function yesorno=cellprescreening_diffusivity(phall,alllabs,Iloci,Ihu,bot,top,left,right);
    % take out brightfield images and stack them vertically
    Ishow2=[];
    for i=1:size(phall,3);
        Ibfcrop = double(phall(bot:top,left:right,i));
        Ibfbg=min(Ibfcrop(:));
        Ibf8bit=Ibfcrop-Ibfbg;
        Ibf8bit=Ibf8bit*255./max(Ibf8bit(:));
        % now take out the Labeled files and crop
        Ilabcrop = double(alllabs(bot:top,left:right,i));
        Ilabbg=min(Ilabcrop(:));
        Ilab8bit=Ilabcrop-Ilabbg;
        Ilab8bit=Ilab8bit*255./max(Ilab8bit(:));
        % now take out the loci files and crop
        Ilocicrop = double(Iloci(bot:top,left:right,i));
        Ilocibg=min(Ilocicrop(:));
        Iloci8bit=Ilocicrop-Ilocibg;
        Iloci8bit=Iloci8bit*255./max(Iloci8bit(:));
        % now take out the hu files and crop
        Ihucrop = double(Ihu(bot:top,left:right,i));
        Ihubg=min(Ihucrop(:));
        Ihu8bit=Ihucrop-Ihubg;
        Ihu8bit=Ihu8bit*255./max(Ihu8bit(:));
    % make stack for display
    Ishow1=cat(1,Ibf8bit,Ilab8bit,Iloci8bit,Ihu8bit);
    Ishow2=cat(2,Ishow2,Ishow1);
    end
    dipshow(Ishow2);
    % dipshow(lab_local,'Labels');
    quest_dualsep = questdlg('Is this a stable single-chromosome cell?','Prescreening','Yes','No','Yes');
    switch quest_dualsep
        case 'Yes'
                close(gcf);
                yesorno=1;

          case 'No'
                 close(gcf);
                 yesorno=0;
    end
end
