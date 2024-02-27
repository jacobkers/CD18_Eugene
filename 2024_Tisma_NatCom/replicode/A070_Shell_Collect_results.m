function A070_Shell_Collect_results
%JWJK_A:-------------------------------------------------------------------
%Shell 

%Description: Shell program to run sequential steps of main code for
%replication analysis. Runs all main programs in set order; former runs can be
%re-used by setting the switches 'if 1'

%input: prior 'crop code' analysis data.

%output:  new .mat database is generated and expanded per step.

%Reference: CD lab, project Sandro, written by Jacob Kers 2018-20
%:JWJK_A-------------------------------------------------------------------

usr='Jacob';
batchrunindices=[0];               %'VersionTest'; 4-chan
% batchrunindices=[1];  %'Raman1';   %5-channel Mukbef with Raman
% batchrunindices=[23]; %initval.expi='Test'; 
batchrunindices=[24.1]; %as 24, test for MukB symmetry,  %same input data as 24, remote paths, add dont data 
batchrunindices=[23 24 25 26 27 28 29 34]; %MukB paper 2022 re-assess data
batchrunindices=[23 24 25 26 27    29 34]; %MukB paper 2022 re-assess data
batchrunindices=[23 24 25 26 27 29 34]; %MukB paper 2022 re-assess data
%batchrunindices=[28]; 

Lb=length(batchrunindices);
dna_results=zeros(Lb,2);
mukbef_results=zeros(Lb,2);
N_okcells=zeros(Lb,1);
patchblock=NaN*zeros(20,1);
for ii=1:Lb
    ii
    bix=batchrunindices(ii);
    disp(['batchrundindex:', num2str(bix)]);
    initval=A000_Repli_Init(bix,usr);     
    [numdat,textdat]=xlsread(strcat(initval.pth_repli,initval.DirSep,initval.expi,'_A062_symmetries.xlsx'),'Summary');
    numdat=[numdat ; patchblock];
    legendz{ii,1}=Replace_underscores(initval.expi);
    dna_results(ii,:)=[numdat(3) numdat(7)]; %median DNA symmetry by area /intensity
    mukbef_results(ii,:)=[numdat(5) numdat(9)]; %median mukbef symmetry  by area /intensity 
    N_okcells(ii)=numdat(1);
end

%plotting:
close all;
figure;
%by area
subplot(2,3,1);
for ii=1:Lb
    bar(ii, dna_results(ii,1)); hold on;    
end
title('Dna by area');
xlabel('exp. index');
ylabel('median of symmetry')
ylim([0 0.7]);

subplot(2,3,2);
for ii=1:Lb
    bar(ii, mukbef_results(ii,1)); hold on;    
end
title('Mukbef by area');
xlabel('exp. index');
ylabel('median of symmetry')
ylim([0 0.7]);

%by intensity
subplot(2,3,4);
for ii=1:Lb
    bar(ii, dna_results(ii,2)); hold on;    
end
title('Dna by intensity');
xlabel('exp. index');
ylabel('median of symmetry')
ylim([0 0.7]);

subplot(2,3,5);
for ii=1:Lb
    bar(ii, mukbef_results(ii,2)); hold on;    
end
title('Mukbef by intensity');
xlabel('exp. index');
ylabel('median of symmetry')
ylim([0 0.7]);

subplot(1,3,3)
for ii=1:Lb
    bar(ii, N_okcells(ii)); hold on;    
end
title('Number of used cells');
legend(legendz, 'Location','NorthOutside');
dum=1;

