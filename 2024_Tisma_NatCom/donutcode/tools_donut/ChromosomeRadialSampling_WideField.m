function Chromosome=ChromosomeRadialSampling_WideField(im,QI, sho);
%JWJK_B:-------------------------------------------------------------------
%Creation of density curves by radial mapping of the chromatin pattern 
%
%Summary: This function analyzes the semi-circular chromosome pattern via radial
%mapping. From this map, 'density curves' are poduced of relative genomic
%content vs. angle. Other properties, such as peak intensity and FWHM width
%are stored as well.
%
%Approach: First it performs COM centering routine on as-is chromosome pattern;
%from this center, radial mapping is performed. Each radial section yields a radial
%profile, which is subjected to peak analysis (position,content, width etc.
%the position of the outer edge of these profiles yoields an 'edge' contour line
%This center-of-mass of this contour line is then used to repeat the
%whole sampling routine. This refinement is performed three times,
%
%Input: image, quadrant interpolation settings, 
%
%Output: structure 'Chromosome' containing density curves and other
%further use
%:JWJK_B-------------------------------------------------------------------
    
[rr,cc]=size(im);

%1) Get COM center of pattern (first estimate)
%QI.maxradius=rr/2;
[xm,ym,~,~,rG]=TrackXY_by_COM_2Dmoment(im);

%2)  Get the radial profiles (running from 0 to 2pi), build a contour line
%and use the coordinates of this line to re-find the COM center
for tt=1:3 
    [rr,cc]=size(im);
     [XX,YY]=meshgrid(1:cc,1:rr);
     Xsamplinggrid=QI.X0samplinggrid+xm;
     Ysamplinggrid=QI.Y0samplinggrid+ym;
     allprofiles=(interp2(XX,YY,im,Xsamplinggrid,Ysamplinggrid,'Linear',0)); 
     %pcolor(flipud(fliplr(im)));colormap bone, shading flat; axis equal; axis tight
     angles=mod(QI.angles,2*pi);
     [ChromosomeAnnularAxis,ix]=sort(angles);
     Chromosome.RadialMap=allprofiles(ix,:); 
     Chromosome.AnnularAxis=ChromosomeAnnularAxis;
     Chromosome.RadialAxis=QI.radii;
      [aa,~]=size(Chromosome.RadialMap);
     Chromosome.MappingIndex=(1:aa);
 
    %3) Add Cell edge and peak detection--------------------------------
     Chromosome=Get_ChromosomeEdges(Chromosome,QI);       
    %4) Get Cartesian Coordiantes
     angles=Chromosome.AnnularAxis;
     polarmax=Chromosome.PolarContourMax;
     polaredge=Chromosome.PolarContourEdge; 
     CX=xm+polaredge.*cos(angles);
     CY=ym+polaredge.*sin(angles);
     MX=xm+polarmax.*cos(angles);
     MY=ym+polarmax.*sin(angles);
     
            
     %update xm, ym
    [xm,ym,theta,ecc]=JKD2_XY_calculate2Dmomentpoints([CX CY],0); 
end
Chromosome.xCOM=xm;
Chromosome.yCOM=ym;
Chromosome.RadGyr=rG;

Chromosome.CartesianContourMax_X=MX;
Chromosome.CartesianContourMax_Y=MY; 
Chromosome.CartesianContourEdge_X=CX;
Chromosome.CartesianContourEdge_Y=CY;

Chromosome.Xsamplinggrid=Xsamplinggrid;
Chromosome.Ysamplinggrid=Ysamplinggrid;

                            
                              
%% get total lengths
    ML=sum(((MX(2:end)-MX(1:end-1)).^2+(MY(2:end)-MY(1:end-1)).^2).^0.5);
    CL=sum(((CX(2:end)-CX(1:end-1)).^2+(CY(2:end)-CY(1:end-1)).^2).^0.5);
    Chromosome.TotalContourLength=CL;
    Chromosome.TotalMaxPeakLength=ML;                            
                              
 if 0
            figure(4);
            MappingIndex=(1:aa);
            VI=Chromosome.ValidIdx;
            cr=QI.radialoversampling
            pcolor(Chromosome.RadialMap); colormap bone; shading flat; hold on;
            plot(cr*Chromosome.PolarContourMax(VI), MappingIndex(VI), 'ro', 'Linewidth', 2);
            plot(cr*Chromosome.PolarContourEdge(VI), MappingIndex(VI), 'bo', 'Linewidth', 1);
            plot(cr*Chromosome.PolarContourRim(VI), MappingIndex(VI), 'yo', 'Linewidth', 1);
            title('radial map of brightfield with bf-edge overlay');
            xlabel('radial sampling steps');
            ylabel('angular sampling steps');
            [~]=ginput(1);
            close(gcf);
        end                             
                              
function Chromosome=Get_ChromosomeEdges(Chromosome,QI);
    %Analyze radial map to find edges etc.      
    RadialMap_ori=Chromosome.RadialMap;    
    [rr,cc]=size(RadialMap_ori);       
    RadialMap=RadialMap_ori;
    Radialweight=1:cc;
    %1) Get some general properties using a 'softened' version
        StetsonMask=BuildStetsonMask(RadialMap); %0.5
    
        MaskedRadialMap=RadialMap_ori.*StetsonMask;
        Chro_Stds=nanstd(MaskedRadialMap'); 
        Chro_Std=nanmedian(std(Chro_Stds));
        GapTreshold=0.02*(4*Chro_Std);
        gapidx=find(Chro_Stds<GapTreshold);  %detect gaps
        nongapidx=find(Chro_Stds>=GapTreshold);  %detect gaps
        Gaps=1.0*Chro_Stds<(GapTreshold)';
        
    if 0
        pcolor(RadialMap); colormap hot; shading flat;
        [~]=ginput(1); close(gcf);
    end
    
    %% 2 Find features of typical radial profile
    for ii=1:rr
        prf=MaskedRadialMap(ii,:);
        
        gap=Gaps(ii);
        Profile=Get_ProfileFeatures(prf,Radialweight,gap);
        FWHM(ii)=Profile.FWHM;
        Content(ii)=Profile.Content; 
        MaxC(ii)=Profile.Max; 
        MaxCVal(ii)=Profile.MaxVal;
        OuterEdge(ii)=Profile.OuterEdge;
        MinC(ii)=Profile.Min;
    end
       
    %padding gaps with neighbouring points OPTIONAL
    %note that measures must be taken to also pad gaps on the edges of the
    %loop. 'MaxC contains the NaN values used to spot gaps
    if 1
        MaxC=Do_Periodic_Padding((1:rr),MaxC,MaxC);
        MaxCVal=Do_Periodic_Padding((1:rr),MaxCVal,MaxC);
        OuterEdge=Do_Periodic_Padding((1:rr),OuterEdge,MaxC);
        MinC=Do_Periodic_Padding((1:rr),MinC,MaxC); 
        
        MaxC=JKD1_PRF_smooth(MaxC',4);  %smooth
        OuterEdge=JKD1_PRF_smooth(OuterEdge',4);  %smooth
        MinC=JKD1_PRF_smooth(MinC',4);  %smooth
    end
    
    Chromosome.PolarContourMax=MaxC/QI.radialoversampling;
    Chromosome.PolarContourEdge=OuterEdge/QI.radialoversampling;;
    Chromosome.PolarContourRim=MinC/QI.radialoversampling;;
    Chromosome.PolarContourContent=Content;
    Chromosome.PolarContourPeakVal=MaxCVal;
    Chromosome.PolarContourFWHM=FWHM/QI.radialoversampling;
    Chromosome.Gaps=Gaps;
    Chromosome.GapIdx=gapidx;
    Chromosome.ValidIdx=nongapidx;
    
 
    if 0
    figure;    
        plot(Chro_Stds); hold on;
        plot(gapidx,Chro_Stds(gapidx),'rx');   
        [~]=ginput(1); close(gcf);
    end
    
    
function StetsonMask=BuildStetsonMask(RadialMap);
%Build a flat-topped edge mask
        Maxdropofffraction=0.9; %max 0.5 for minima
        Mindropofffraction=0.1;
        [rr,cc]=size(RadialMap);
        
        %build a global profile
        EstDensityCurve=sum(RadialMap');
        EstDensityCurve=EstDensityCurve/nanmax(EstDensityCurve);
        InvEstDensityCurve=1-EstDensityCurve;
        
        dropofffractioncurve=(Maxdropofffraction-Mindropofffraction)*...
                             InvEstDensityCurve+...
                             Mindropofffraction;
        
        StetsonMask=zeros(rr,cc);
        %fill it with amplitude dependent masking
        for jj=1:rr
            dropofffractionIn=dropofffractioncurve(jj);
            dropofffractionOut=nanmin(dropofffractioncurve);
            dropofflengthIn=ceil(cc*dropofffractionIn);
            dropofflengthOut=ceil(cc*dropofffractionOut/2);
            
            flatlength=cc-dropofflengthIn-dropofflengthOut;
            EdgesIn=(hann(2*dropofflengthIn))';
            EdgesOut=(hann(2*dropofflengthOut))';
            StetsonCurve=[EdgesIn(1:dropofflengthIn) ...
                         ones(1,flatlength)...
                         EdgesOut(dropofflengthOut+1:end)];    
            StetsonMask(jj,:)=StetsonCurve;    
        end
        dum=1;

function Profile=Get_ProfileFeatures(prf,Radialweight,gap);
    %JWJK_C:
    %-------------------------------------------------------------------------
    %This function obtains properties of one radial profile (typically, of 256)
    %first, profile is multiplied with hanning, initial max is determined there. 
    %(to eliminate artefacts from neighbouring cells, especially near 
    %low-intensity ter gaps etc. Then max-search is performed 
    %on the original profile from this point on. A FWHM is obtained around this
    %mask. Results are corrected for the radial oversampling we used and are in
    %image pixel units. As a measure for local intensity, we take the sum of
    %the profile but weighed per radial position to acount for the different
    %number of contributing pixels (outer image pixels should count more). The
    %outer edge is defined as the largest radial position with an intensity above
    %35% of the peak value.
    %------------------------------------------------% Jacob Kerssemakers 2016 
    %:JWJK_C


%Get typical radial features of profile
prf=prf-min(prf);
cc=length(prf); 
if ~gap   
    %1) find first profile max from center pos
    prf_han=prf.*(hann(cc)');
    [~,MaskCenter]=max(prf_han);
    MaskWidth=50;
    %prf_smz=JKD1_PRF_smooth(prf',1);     
    maxidx=Get_GuidedMax(prf,MaskCenter,MaskWidth);
    peakval=prf(maxidx);
    if 0
        figure;    
        plot(prf,'b-'); hold on;
        %plot(prf_smz,'r-');   
        plot(maxidx,peakval,'ko');  
        [~]=ginput(1); close(gcf);
    end
    
    outprf=prf(maxidx:end);  %profile on inside 
    inprf=prf(1:maxidx-1);   %profile on outside
      
    %dark rim-------------------------------------
    [~,ix]=min(outprf);
    if isempty(maxidx)
        dum=1;
    end

    outminidx=ix+maxidx-1;    
    %build edge----------------------------
    sel2=find((prf>0.35*peakval)&(~isnan(prf)));
    edgeidx=max(sel2);
    
    %content, note that we correct for radial weight
    sel=find(prf>0);
    cnt=nansum(prf(sel).*Radialweight(sel));
    
    %width
    sel2=find((prf>0.5*peakval)&(~isnan(prf)));
    FWHM=length(sel2);
    
    Profile.FWHM=FWHM;
    Profile.Content=cnt; 
    Profile.Max=maxidx;
    Profile.MaxVal=peakval;
    Profile.Min=outminidx;       
    Profile.OuterEdge=edgeidx;

else
    Profile.FWHM=NaN;
    Profile.Content=0;
    Profile.MaxVal=0;
    Profile.Max=NaN;
    Profile.Min=NaN;
    Profile.OuterEdge=NaN;    
end 
    


function outdata=JKD1_PRF_smooth(data,span)
%running window 'span' smoothing average; 0=no smoothing.
%perioc boundaries
%JacobKers 2013

if span>0
    hs=ceil(span/2);
    le=length(data);
    outdata=zeros(le,1);
    midpt=ceil(le/2);
    data2=[data(midpt:end) ; data ; data(1:midpt)];
    lex=length(data(midpt:end));
    for i=1:le
    idx=i+lex;
    outdata(i)=mean(data2(idx-hs:idx+hs));
    end
else
    outdata=data; %no smoothing
end

function data_expanded=Expand_elliptoidrange(data);
    %Expand for later interpolation
    [rr,~]=size(data);
    halfidx=ceil(rr/2);
    half1=data(1:halfidx-1,:);
    half2=data(halfidx:end,:);
    half1(:,1)=half1(:,1)+2*pi;
    half2(:,1)=half2(:,1)-2*pi;
    data_expanded=[half2 ; data ; half1];



           
 
       



