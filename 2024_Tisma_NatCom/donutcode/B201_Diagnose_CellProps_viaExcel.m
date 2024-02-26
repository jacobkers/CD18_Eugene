function B201_Diagnose_CellProps_viaExcel
close all;
batchrunindex=-101.1;
initval=A000__WF_Get_JacobPathsandExperiments(batchrunindex);


analyze_cellshapes=0;
analyze_clusters=1;
channel='c4';

%% 2 cell shape analysis (try to mimick filter props)
if analyze_cellshapes
    %read excel. the 'numdat' rows are shifted one, column indexes are identical
    %for re-writing, use first numdat row to detrmine text or numeric input
    %(for text, numeric is empty)
    [numdat,txtdat]=xlsread(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A010_Cell_BasicGeometries.xlsx'));   
    HeadersIN=txtdat(1,:); 
    ValuesIN_num=numdat;
    HeadersINPerCell=txtdat(1,:); 
    ValuesINPerCell=numdat(3:end,:);
    
    % {'index'                           }
    % {'label'                           }
    % {'chromosome_peak_length'          }
    % {'chromosome_outercontour_length'  }
    % {'chrom_peak_averageradius'        }
    % {'chrom_peak_stdradius'            }
    % {'chrom_peak_maxradius'            }
    % {'chrom_peak_minradius'            }
    % {'chrom_outercontour_averageradius'}
    % {'chrom_outercontour_stdradius'    }
    % {'chrom_outercontour_maxradius'    }
    % {'chrom_outercontour_minradius'    }
    % {'chrom_radiusofgyration'          }
    % {'chrom_averageFWHM'               }
    % {'chrom_stdFWHM'                   }
    % {'chrom_maxFWHM'                   }
    % {'chrom_minFWHM'                   }
    % {'cell wall length'                }
    % {'cellwall_averageradius'          }
    % {'cellwall_stdradius'              }
    % {'cellwall_maxradius'              }
    % {'cellwall_minradius'              }
    % {'rawcounts_chro'                  }
    % {'rawcounts_rfp'                   }
    % {'rawcounts_cfp'                   }
        
    %columns
    cell_wall_length= ValuesIN_num(1:end,find(strcmp(HeadersIN,'cell wall length')));
    cellwall_averageradius= ValuesIN_num(1:end,find(strcmp(HeadersIN,'cellwall_averageradius')));
    cellwall_stdradius= ValuesIN_num(1:end,find(strcmp(HeadersIN,'cellwall_stdradius')));
    cellwall_maxradius= ValuesIN_num(1:end,find(strcmp(HeadersIN,'cellwall_maxradius')));
    cellwall_minradius= ValuesIN_num(1:end,find(strcmp(HeadersIN,'cellwall_minradius')));    
    %derived
    roundness=cell_wall_length./cellwall_averageradius/(2*pi);
    elongation=(cellwall_maxradius-cellwall_minradius)./cellwall_averageradius;
    diameter=2*cellwall_averageradius;
    
    %roundness
    axz_roundness=linspace(0,4,20);
    hist_roundness=hist(roundness,axz_roundness);
    subplot(2,3,1);
        bar(axz_roundness,hist_roundness,'r');
        xlabel('perim/2piR')
        ylabel('counts');
        title('roundness');
        legend(Replace_underscores(initval.expi));
    subplot(2,3,4);
        bar(axz_roundness,cumsum(hist_roundness),'r');
        xlabel('perim/2piR')
        ylabel('total counts');
        title('roundness');
    %elongation    
    axz_elongation=linspace(0,4,20);
    hist_elongation=hist(elongation,axz_elongation);
    subplot(2,3,2);
        bar(axz_elongation,hist_elongation,'b');
        xlabel('(Rmax-Rmin)/Rav')
        ylabel('counts');
        title('elongation');
    subplot(2,3,5);
        bar(axz_elongation,cumsum(hist_elongation),'b');
        xlabel('(Rmax-Rmin)/Rav')
        ylabel('total counts');
        title('elongation');
    %diameter
    axz_diameter=linspace(0,50,20);
    hist_diameter=hist(diameter,axz_diameter);
    subplot(2,3,3);
        bar(axz_diameter,hist_diameter,'m');
        xlabel('av.diameter')
        ylabel('counts');
        title('diameter');
    subplot(2,3,6);
        bar(axz_diameter,cumsum(hist_diameter),'m');
        xlabel('av.diameter')
        ylabel('total counts');
        title('diameter');
end

dum=1;



%% 2 cell cluster analysis
if analyze_clusters
    %read excel. the 'numdat' rows are shifted one, column indexes are identical
    %for re-writing, use first numdat row to detrmine text or numeric input
    %(for text, numeric is empty)
    [numdat,txtdat]=xlsread(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A050_ClusterColumnReport_',channel,'.xlsx'));
    [numdat2,txtdat2]=xlsread(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A050_Cell_ClusterReport_',channel,'.xlsx'));   
    HeadersIN=txtdat(1,:); 
    ValuesIN_num=numdat;
    HeadersINPerCell_1=txtdat2(1,:); 
    HeadersINPerCell_2=txtdat2(2,:); 
    ValuesINPerCell=numdat2(3:end,:);

    %somespecific column actions
    CellIndex=ValuesIN_num(1:end,find(strcmp(HeadersIN,'index')));
    Contents= initval.genomelength*ValuesIN_num(1:end,find(strcmp(HeadersIN,'content')));
    Areas=(initval.nmperpixel)^2*ValuesIN_num(1:end,find(strcmp(HeadersIN,'area')));
    RadGyr3=initval.nmperpixel*ValuesIN_num(1:end,find(strcmp(HeadersIN,'diameter_ContourBW_Equiv')));
    ClusterNoPerCell=ValuesINPerCell(1:end,find(strcmp(HeadersINPerCell_2,'cluster number')));

    MxA=nanmax(Areas);
    axz=linspace(0,MxA,10);
    HistArea=hist(Areas,axz);
    MxC=nanmax(Contents);
    AxzC=linspace(0,MxC,10);
    HistContent=hist(Contents,AxzC);
    MxN=nanmax(ClusterNoPerCell);
    AxzN=0:1:MxN;
    HistClusterNo=hist(ClusterNoPerCell,AxzN);

    subplot(1,3,1);
        bar(AxzN,HistClusterNo,'r');
        xlabel('nluster number per cell)')
        ylabel('counts');
        xlim([1 8]);
        ylim([0 135]);
    subplot(1,3,2);
        bar(AxzC,HistContent,'b');
        xlabel('cluster DNA content(bp)')
        ylabel('counts');
        xlim([-1000 4000]);
    subplot(1,3,3);
        plot(Contents, RadGyr3,'ko', 'MarkerSize',5, 'MarkerFaceColor','y');
        xlabel('cluster DNA content(bp)')
        ylabel('equivalent diameter, nm');
        ylim([0 1000]);
end







