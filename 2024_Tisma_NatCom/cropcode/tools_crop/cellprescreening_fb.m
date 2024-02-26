function  [addlab,labB,labreplace]=cellprescreening_fb(labA,IchA,labind,addlab,MagDisplay,autorunit);
rounds=0;
labB=[];
    while rounds==0
        if isempty(labB)==0; labA=labB;
        end
        Ishow2=[];
    Ilab8bit=labA*255./max(labA(:));
    Ilab2bit=labA>0;
    Ilab2bitdb=double(Ilab2bit);
     for i=1:size(IchA,3) % scale all images and turn into 8 bit
         IchB=IchA(:,:,i);
         Ibg=min(IchB(:));
         I8bit=IchB-Ibg;
         I8bit=I8bit*255./max(I8bit(:));
         if i==1
             I8bitA=I8bit;
             I8bitA(Ilab2bit==1)=80;
             Ishow1=[I8bitA Ilab8bit];  
         else
         I8bitB=I8bit.*Ilab2bitdb;
             Ishow1=[I8bit I8bitB];
         end
         Ishow2=cat(1,Ishow2,Ishow1);
     end
     Ishow3=imresize(Ishow2,MagDisplay,'bilinear');
     
        if ~autorunit
        dipshow(Ishow3);
        quest_step1 = questdlg('Hey dude, what do you want?','Prescreening','Just fine','Wanna fix it','Trash it','Just fine');
            switch quest_step1
                case 'Just fine'
                    labreplace=2; labB=labA; rounds=1;
                    close(gcf);
                case 'Trash it'
                    labreplace=1; labB=labA<0;rounds=1;
                    close(gcf);
                case 'Wanna fix it'
                    labreplace=2;rounds=0;
                    quest_step2 = questdlg('So Shall we?','fix option','Patch it','Slice it','Move on','Patch it');
                    switch quest_step2
                       case 'Patch it'
                           boundtype=1;
                           [addlab,labB] = add_boundaryfb(Ilab2bit,boundtype,labind,addlab,MagDisplay);
                           close(gcf);
                        case 'Slice it'
                            boundtype=0;
                           [addlab,labB] = add_boundaryfb(Ilab2bit,boundtype,labind,addlab,MagDisplay);
                           close(gcf);
                        case 'Move on'
                            labB=labA;rounds=1;
                            close(gcf);
                    end
            end
        else  %pass all cells
           labreplace=2; labB=labA; rounds=1;
           pause(0.1);
            close(gcf); 
        end
    end
end
