function [Stack,workpic]=Get_StackMovieFromDefocusSeries(initval,cellno,channelname,WorkPicOption);
    %build stacks from z-focus planes  as in
    %'20160126bn2179a22t01rfp_cell1z1.tif';
    switch initval.storagemodus
        case 1,             
            datapath=strcat(initval.allchannelpath,channelname,'\'); 
            cell_label=strcat('*cell',cellno);
             cd(datapath); 
             AllFiles=dir(strcat('*',cell_label,'*.tif'));
             filelabel='20160126bn2179a22t*';
        case 2 
        datapath=strcat(initval.allchannelpath,'cell',cellno,'\'); 
            switch channelname
                case 'bf', cell_label=strcat('*c4cell',cellno);
                case 'cfp', cell_label=strcat('*c1cell',cellno);
                case 'rfp', cell_label=strcat('*c2cell',cellno);
                case 'yfp', cell_label=strcat('*c3cell',cellno);
                case 'bf', cell_label=strcat('*c4cell',cellno);
            end
            cd(datapath); AllFiles=dir(strcat(cell_label,'.tif'));
            filelabel='bn2179_a22_xy*_t*';
    end
 
    
    %1) get all planes in 4D matrix
    WhatPlane='FocalPlane';
    AllPlanes=3;
    AllTimes=11;
    
   
   
    A0 = imread(strcat(datapath,AllFiles(1).name));
   [rr,cc]=size(A0);
    SuperStack=zeros(rr,cc,AllTimes,AllPlanes);
    for planeno=1:AllPlanes
            switch initval.storagemodus
                case 1, z_label=strcat(filelabel,channelname,cell_label, 'z',num2str(planeno),'.tif');
                case 2, z_label=strcat(filelabel,'z',num2str(planeno), cell_label,'.tif');
            end
            %this loads the z-planes of one frametime    
            Files=dir(z_label);            
            for kk = 1:AllTimes
                A = double(imread(strcat(datapath,Files(kk).name)));
                SuperStack(:,:,kk,planeno)=A;
               % [~]=ginput(1); pause(0.1);
            end  
    end
    cd(initval.mainpath);  
    %get main plane  (assuming largely planar features)
    
    
    %2) per time; get ONE plane
    Stack=zeros(rr,cc,AllTimes);
    for tt = 1:AllTimes
        ZStackT=zeros(rr,cc,AllPlanes) ;
        for pl=1:AllPlanes
            ZStackT(:,:,pl)=squeeze(SuperStack(:,:,tt,pl));
        end        
        switch WhatPlane
            case 'MeanProject',
                BestPlane=squeeze(mean(ZStackT,3));
            case 'FocalPlane', 
                %get main plane  (assuming largely planar features)
                stcurve=squeeze(std(std(ZStackT)));
                [~,MainPlane]=max(stcurve);
                BestPlane=double(ZStackT(:,:,MainPlane));  
        end
        Stack(:,:,tt)=BestPlane;
    end
    
  
    %3 get the 'workpicture'
    
        switch WorkPicOption
            case 'MeanProject',
                workpic=squeeze(mean(Stack,3));
            case 'FocalPlane', 
                %get main plane  (assuming largely planar features)
                stcurve=squeeze(std(std(Stack)));
                [~,MainPlane]=max(stcurve);
                workpic=double(Stack(:,:,MainPlane));  
            case 'First', 
                workpic=double(Stack(:,:,1));
        end

    %------------------------------------------------------------    
