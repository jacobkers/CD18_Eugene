function [Stack,workpic]=Get_Stack(initval,cellno,channelname,WorkPicOption);
    
    
    switch initval.storagemodus
        case 1,             
            datapath=strcat(initval.allchannelpath,channelname,'\'); 
            labl=strcat('*cell',cellno);
            cd(datapath); Files=dir(strcat(labl,'.tif'));
        case 2 
            datapath=strcat(initval.allchannelpath,'cell',cellno,'\'); 
            switch channelname
                case 'bf', labl=strcat('*C4cell',cellno);
                case 'yfp', labl=strcat('*C3cell',cellno);
            end
                cd(datapath); Files=dir(strcat(labl,'.tif'));
                dum=1;
    end
    
    
    cd(initval.mainpath);  
    num_images=length(Files);
    A0 = imread(strcat(datapath,Files(1).name));
    [rr,cc]=size(A0);
    Stack=zeros(rr,cc,num_images);
    for kk = 1:num_images
        A = double(imread(strcat(datapath,Files(kk).name)));
        Stack(:,:,kk)=A;
       % [~]=ginput(1); pause(0.1);
    end  

    %get main plane  (assuming largely planar features)
    
    %inpic=double(Stack(:,:,1));    
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
