function A401_Correlation_Plotter
%JWJK_A:-------------------------------------------------------------------
%Title: matrix plotting of results
%
%Summary: This function perfoms matrix plotting of  data processed by A400
%
%Input: mat files in specified output dirspecified
%Output: matrix plots

%Reference: Cees Dekker Lab, Project: MinED; researcher Grzegorz Pawlik; 
%code designed & written Jacob Kerssemakers 2016 
%:JWJK_A-------------------------------------------------------------------

clear all; close all;
addpath(genpath(pwd));

outdirname=swap_path('D:\jkerssemakers\Dropbox\CD_Data_out\2016_Greg\20180903 Correlation\');
experimentseriesname='movies_4';  %the name used for this series




load(strcat(outdirname,experimentseriesname,'_matlaboutput\A400_Correlation_results_tabulated.mat'),...
    'heights', 'EDratios', ...
    'xtcor_decaylength_abs_matrix', ...
    'xtcor_decayperiodicityfactor_matrix', ...
    'xtcor_PK1radii_matrix',...
    'xtcor_domainL_matrix', ...
    'xycor_decaylength_abs_matrix', ...
    'xycor_decayperiodicityfactor_matrix', ...
    'xycor_PK1radii_matrix');

    
figure(2);

plotits.panelidx=1;
plotits.title='wavelength';
plotits.xlabel='height index';
plotits.ylabel='ED-ratio index';
heatmapit(EDratios,heights,xycor_PK1radii_matrix,plotits); 


function heatmapit(EDratios,heights,data_of_interest,plotits)                        
%build graphs of corelation peaks as function of height, per ED ratio;
subplot(2,2,plotits.panelidx); % channel

pcolor(data_of_interest); shading flat; colormap hot
title(plotits.title);
xlabel(plotits.xlabel);
ylabel(plotits.ylabel);

