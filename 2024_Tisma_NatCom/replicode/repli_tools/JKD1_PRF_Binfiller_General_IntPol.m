function binresults=JKD1_PRF_Binfiller_General_IntPol(Xdata2D,Ydata2D,binax);
% This function allocates data to a series of bins and averages the result;

% If the number of points in the X-axis is small (< bins) and normalized, direct filling is irregular. 
% In this version, the data is first interpolated to ensure each bin is
% filled once and only once per trace
% Use this version if you are  not interested in counts, but in averages


%input: X data, Y data: 2D array of index, values. Xdata en Y data is assumed to be
%concatenated; X data is assumed to increase in time per sub-trace
%output:
    % binresults.losig=lo;          %one sigma below the average, per bin
    % binresults.av=av;             % the average, per bin
    % binresults.hisig=hi;          %one sigma above the average, per bin
    % binresults.binaxis=binaxnw;   %values of the non-empty bins, per bin
    % binresults.binaxis_all=binax2D(:); %value of  bin, per point
    % binresults.binaxis_all_scat=scatax2D(:); %random number to spread bin points (for plotting purposes)
    % binresults.bincounts=bincounts;
    
% JacobKers 2013
%------------------------------------------------------


scatwidth=0.7;

%---------------------------------------------------
if nargin<3
    close all
    minbin=-0.1;
    maxbin=1.1;
    bins=20;
    binax=linspace(minbin,maxbin,bins);
    binsize=binax(2)-binax(1);
    %fake data
    [Xdata2D,Ydata2D]=CreateFakeData;
end
%----------------------------------------------------
bins=length(binax);
binsize=binax(2)-binax(1);
minbin=binax(1);
maxbin=binax(end);

Idxdata=Xdata2D(:,1);
Xdata=Xdata2D(:,2);
Ydata=Ydata2D(:,2);

%remove outliers looking at the whole dataset
[flag,~]=JKD1_PRF_outlier_flag(Ydata,3,0.8,'all',0);
sel=find(flag==1);
Xdata=Xdata(sel);
Ydata=Ydata(sel);
Idxdata=Idxdata(sel);
% -------------------------------

%First, measure the number of separate datasets by looking at 'negative time
%jumps'; define bincollector matrix
dif=Idxdata(2:end)-Idxdata(1:end-1);
sel=(find(dif<0)); 
traceno=length(sel)+1;
tracestart=([0 sel' length(Idxdata)])';
bins_L=length(binax);
bincollector=zeros(bins_L,traceno)*NaN;  %define bins

%Next, transform each dataset and fill bins

for i=1:traceno    
    x=Xdata(tracestart(i)+1:tracestart(i+1))';
    y=Ydata(tracestart(i)+1:tracestart(i+1))';   
    [xq,ixq,~]=unique(x);   %remove double entries from earlier steps
    yq=y(ixq);
    xip=binax';
    yip=interp1(xq,yq,xip,'linear', NaN);
    bincollector(:,i)=yip;
    if 0; 
    plot(x,y,'o-'); hold on;
    plot(xip,yip,'ro'); hold off;    
    [~]=ginput(1);
    end 
    dum=1;
end


%Define a 'scatter' axis and an axis for each individual entry in the bins
scatax2D= scatwidth*binsize*(rand(bins_L,traceno)-0.5);
binax2D=repmat(binax',1,traceno);

%clean the peak data from outliers and get averages-----------------
av=zeros(bins_L,1);
sig=zeros(bins_L,1);

for i=1:bins_L    
sel=~isnan(bincollector(i,:));
onebindata=bincollector(i,sel);
    if length(onebindata)>0
    %Detect 'inliers' and refine selection--------------------
    [flag,cleandata]=JKD1_PRF_outlier_flag(onebindata,2,0.8,'positive',0);
    %data,tolerance,sigchange,how,sho, binz
    sel2=find(flag==1);
    nonsel=find(flag==0);
    bincollector(i,sel(nonsel))=NaN;
    
    %Get averages and standard deviation of final selection
    av(i)=nanmean(onebindata(sel2));
    sig(i)=nanstd(onebindata(sel2));
    else
    av(i)=0;
    sig(i)=0;
    end
end
%Finally, upper and lower bounds (one sigma)-------------------------------
lo=av-1*sig;
hi=av+1*sig;

binresults.losig=lo;          %one sigma below the average, per bin
binresults.av=av;             % the average, per bin
binresults.hisig=hi;          %one sigma above the average, per bin
binresults.binaxis=binax;   %values of the non-empty bins, per bin
binresults.binaxis_all=binax2D(:); %value of  bin, per point
binresults.binaxis_all_scat=binax2D(:)+scatax2D(:); %random number to spread bin points (for plotting purposes)
binresults.binvalues_all=bincollector(:);
binresults.bincounts=traceno;


%Plot menu--------------------------------------------------------------
%1a)  total&peaks per bacterium per time point
if nargin<4
figure;
plot(binax2D+scatax2D,bincollector, 'bo', 'MarkerSize', 1.5); hold on;
title(strcat('Binned intensities-'));
xlabel('frames');
ylabel('Intensity per length');
axis([minbin max(binax) 0 max(Ydata)]);

%1b totals peaks av-lo-hi
plot(binax,av,'--bs','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','w',...
                'MarkerSize',5); hold on;
plot(binax,lo,'bo-','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','w',...
                'MarkerSize',5); hold on;
plot(binax,hi,'bo-','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','w',...
                'MarkerSize',5); hold on; 
%--------------------------------------------------------------------------
end

function [Xdata,Ydata]=CreateFakeData;
%This function generates traces with length varying randomly around an average
%input: mininimum value X, maximum value X.
%output: Ntrace of these concatenated, X and Y
%It shows the '0' and '1' gaps artefacts
%JacobKers 2013

Ntrace=200;
Avpointspertrace=20;
lox=-0.3;
hix=1.8;
loy=20;
hiy=40;
Xdata=[]; Ydata=[];

for i=1:Ntrace
tracelength=ceil(Avpointspertrace*(1+0.5*(rand(1,1)-0.5))); %50% plus min variation
x=linspace(lox,hix,tracelength)';
y=linspace(loy,hiy,tracelength)'+0.5*(hiy-loy)*(rand(1,tracelength)-0.5)';
if 0; 
    plot(x,y,'o-');
    %[~]=ginput(1);
end
idx=(1:tracelength)';
Xdata=[Xdata ; [idx x]]; 
Ydata=[Ydata ; [idx y]];
end
