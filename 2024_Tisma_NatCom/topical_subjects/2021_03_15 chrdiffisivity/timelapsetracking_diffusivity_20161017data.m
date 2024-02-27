%% This code takes whatever time-lapse stacks and measure the positions of the cells
%% chromosomes and foci. Data are names *c*.tif. 
%% Note that this code restrains only one focus per foci channel

dircrop='/Users/fabaiwu/Documents/work/drafts/ChrCompactionSegregation/Figure2/material/diffusivity/20161017shorterm/allcropped';

cd(dircrop);
dirc1=dir('*c1.tif');
chrch='c2'; % leave it to '' if none
focich=cat(1,'c3','c4'); % leave it to {} if none
numfoci=size(focich,1);
% cellstrs={};
for i=110:length(dirc1)
    if mod(i,25)==0; 
        i 
    end
    cellstr=dirc1(i).name(1:21);
    celldatname=[cellstr '_data.mat'];
    % cellstrs={cellstrs,cellstr};
    if numel(chrch)>0
       Ichrname=[cellstr '_' chrch '.tif'];
    end
    if numfoci>0
        Ifocinames='';
        for j=1:numfoci
            fociname1=cat(2,cellstr,'_',focich(j,:),'.tif');
            Ifocinames=cat(1,Ifocinames,fociname1);
        end
    end
    Iinfo=imfinfo(dirc1(i).name);
    tfrs=length(Iinfo);
    celldata=[];chrdata=[];focidata=[];
    for t=1:tfrs
        phase1=double(imread(dirc1(i).name,t));
        % [imh,imw]=size(phase1);
        Iph = double(gaussf(phase1,10,'best') - gaussf(phase1,2,'best'));% Using Gaussian to smooth the image analyzed
        % thresholding to find cells
        thres=0.5*max(Iph(:));
        Iph=Iph>thres;
%        Iph=Iph-bpropagation(newim(Iph,'bin'),Iph,0,2,1);
        Iph=brmedgeobjs(Iph,2); % get rid of cells cut off at the edge.
        Iph1=imfill(logical(Iph),'holes');
        lab1=label(Iph1,2,300,10000); % set threshold for the cell size        
        labph=double(lab1);
        modeph=mode(labph(labph>0));
        if isnan(modeph)==1
            continue
        else
            Iph2=labph==mode(labph(labph>0)); % find the biggest object

            data=measure(Iph2,[],{'size','feret','center','minimum','maximum'}); % minimum/maximum are shown in [x;y]
            % ferets data: maxFeret, minFeret, perpenFeret, maxFangle, minFangle (note, angles are -2pi to 2pi)
            celldat1=cat(2,t-1,double(data.size)',double(data.center)',double(data.minimum)',double(data.maximum)',double(data.feret)');
            %% Now open the HU channel
            if numel(chrch)>0
                chr1=double(imread(Ichrname,t));
                chr1=chr1-mode(chr1(:)); % background subtraction
                chr1(chr1<0)=0;
                chr2=gaussf(chr1,2,'best') - gaussf(chr1,10,'best'); 
                thres=0.3*max(chr2(:));
                chr2=chr2>thres;
                lab2=label(chr2,2,8,2000);
                lab2=lab2>0;
                labchr=double(lab2).*labph;
                [y2,x2]=find(labchr==1);
                Ichr=chr1.*labchr; % 
                intchr=[0 0];
                for k=1:numel(y2)
                    intchr=intchr+[x2(k) y2(k)]*Ichr(y2(k),x2(k));
                end
                sumchr=sum(Ichr(:));
                areachr=sum(labchr(:));
                meanchr=sumchr/numel(y2); % mean chromosome intensity
                centchr=intchr./sumchr; % center of mass
                chrdat1=cat(2,t-1,areachr,centchr,min(x2),min(y2),max(x2),max(y2),meanchr);
            end
            % now work on foci
            if numfoci>0
                focidat1=[];
                for ff=1:numfoci
                    focus1=double(imread(Ifocinames(ff,:),t));
                    focus1=focus1-mode(focus1(:)); % background subtraction
                    focus1(focus1<0)=0;
                    focus2=double(gaussf(focus1,2,'best') - gaussf(focus1,10,'best')); 
                    maxV=max(focus2(:));
                    [Ym,Xm]=find(focus2==maxV); % find the location of the brightest spot
                    thres=0.3*maxV;
                    focus2=focus2>thres;
                    lab3=label(focus2,2,4,500);
                    labfocus=double(lab3).*labph; 
                    labfocus=labfocus==mode(labfocus(Ym(1),Xm(1))); % the spot with the brightest pixel
                    [y2,x2]=find(labfocus==1);
                    Ifocus=focus1.*labfocus; % 
                    intfocus=[0 0];
                    for k=1:numel(y2)
                        intfocus=intfocus+[x2(k) y2(k)]*Ifocus(y2(k),x2(k));
                    end
                    sumfocus=sum(Ifocus(:));
                    areafocus=sum(labfocus(:));
                    meanfocus=sumfocus/numel(y2); % mean focus intensity
                    centfocus=intfocus./sumfocus; % center of mass
                    focidat1=cat(3,focidat1,cat(2,t-1,areafocus,centfocus,min(x2),min(y2),max(x2),max(y2),meanfocus));
                end
            end
            celldata=cat(1,celldata,celldat1);
            chrdata=cat(1,chrdata,chrdat1);
            focidata=cat(1,focidata,focidat1); 
        end
    end
    save(celldatname,'celldata','chrdata','focidata');
end
% save('allcellnames','cellstrs');
    
            
        
        