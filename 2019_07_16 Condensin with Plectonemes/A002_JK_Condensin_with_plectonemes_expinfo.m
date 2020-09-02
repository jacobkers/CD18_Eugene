function expinfo=A002_JK_Condensin_with_plectonemes_expinfo(expi,roi);
%JWJK_A:----[add ABCorC*----------------------------------------------------
%Title: %Settings per experiment. 
%Summary: Use imageJ to get these numbers settings for all
%files. the number ''expinfo.channelshift'' is to be gotten from the
%apparent position shift in condnesin and plectoneme kymographs
%Input: index of experiment
%Output: 
%References: CD lab, project Eugene Kim, written by Jacob Kers, 2019
%:JWJK_A-----[add ABCorC*---------------------------------------------------


%expinfo.channelshift=0.6022; %manually per ROI....  
switch expi
    case -1,  expinfo=get_roiproperties_of_2020_08_18_Simulation(roi);
    case 0,   expinfo=get_roiproperties_of_2019_09_02_NumberOfCondensinPerDNA(roi);
    case 1,   expinfo=get_roiproperties_of_2019_07_15_condensin_supercoil(roi);          
    case 2,   expinfo=get_roiproperties_of_2019_07_26_condensin_supercoil_no_ATP(roi);         
    case 3,   expinfo=get_roiproperties_of_2019_09_02_NumberOfCondensinPerDNA(roi);  
    case 4,   expinfo=get_roiproperties_of_2020_01_13_MukBEF_msd_wATP(roi); 
    case 5,   expinfo=get_roiproperties_of_2020_05_05_data_sc_cnd(roi);
    case 11,  expinfo=get_roiproperties_of_Atto_condensin_42kb_nicking(roi); 
    case 12,  expinfo=get_roiproperties_of_WT_condensin_42kb_non_nicking(roi); 
end

function expinfo=get_roiproperties_of_2020_08_18_Simulation(roi); 
        %%simulations, just filling the fields
        expinfo.labelname='Condensin'; %no condensin channel
        expinfo.tres_pk_DNA=0.4;      %used for final peak selection; fraction of first 90% %change this for more conservative peak detection
        expinfo.tres_pk_Cnd=4;       %used for final peak selection; sigmas over noise
        roiprops=[...  %no  x1  y1  x2   y2     wd  drx dry chanshift
                        1   NaN  NaN   NaN   NaN    NaN  NaN  NaN  0; %shift upwards (-) ;leftwards (-)
         ];
         expinfo.endpoints_xy=[roiprops(roi,2:3) ; roiprops(roi,4:5)]; %in image x1 y1 x2 y2   
         expinfo.driftxy=roiprops(roi,7:8); %deltax, deltay %between startframe and endframe
         expinfo.kymowidth=roiprops(roi,6); %adjust with neighbours nearby
         expinfo.channelshift=roiprops(roi,9); %manually per ROI.... 

function expinfo=get_roiproperties_of_WT_condensin_42kb_non_nicking(roi); 
        %'Atto_condensin_22kb_non_nicking\'
        expinfo.labelname='DNA'; %no condensin channel
        expinfo.tres_pk_DNA=0.4;      %used for final peak selection; fraction of first 90% %change this for more conservative peak detection
        expinfo.tres_pk_Cnd=4;       %used for final peak selection; sigmas over noise
        roiprops=[...  %no  x1  y1  x2   y2     wd  drx dry chanshift padit           
               1    16   3  25   119   14   0  0   0 0; 
               2    34   8   8   58    17   0  0   0 0;
               3    29  12  26   82   20    0  0   0 0;
               4   16  10  26  83    15     0  0   0 1;
               5   34  10  14  65   18      0  0   0 0;
               6   10  2   8   65   15      0   0  0 1;
               7   10  4   18  110   11     0   0  0 1;
               8   9   5   25  90    17     0   0 0 0;
               9   21  18  34  101   17     0   0 0 0;
               10  15  15  27  101   17     0   0  0 0;
               11  27  13   16   116  17    0   0  0 0;
               12  10  5   20   88   17     0   0  0 1;
               13  36  9   17   110   17    0   0  0 0;
               14  20  8   26   113   20    0   0  0 1;
               15  31  2   16   50   15     0   0  0 1;
               16  18  12   31  53   20     0   0  0 1;
               17  20  10   16  80   20     0   0  0 0;
               18  24  9   15   113  11    -1   0  0 0;
               19  25  9   9   124  11      -3   0  0 1;
               20  8   9   33   65   14     0   0  0 1;
               21  12  11  34   88   15     0   0  0 1;
               22  16  9  18  88   10       0   0 0 0;
               23  15  8   35   90  15      0   0  0 0;
               24  15  8   17   74   15     0   0  0 0;
               25  12  7   15   66   15     0   0  0 0;
               26  11  9   14   69   17     0   0  0 0;
               27  20  4   20   77   15   -4   0  0 0;
               28  14  2   22  64   20      0   0  0 0;
               29  26  11   14  73   17     0   0  0 0;
               30  27  9   20   87  20      0   0  0 0;
               31  12  21   56   56  20     0   0  0 0;
               32  29  12  61   65   41   -2  -2  0 0;
               33  15  12  24   114   20   0   0  0 0;
               34  13  6   30   85   20     0   0  0 0;
               35  14  12   22   114   20   0   0  0  0;
               36  14  6   19   103   15   0   0  0 0;
               37  14  5   19   104   15   0   0  0 0;
               38  34  4   24  80   10   0   0  0 0;
               39  30  10   31  142   10   0   0  0 0;
               40  12  32   73   115  15   0   0  0 0;
               41  16  6   25   107  10  0   0  0 0;
];
         expinfo.endpoints_xy=[roiprops(roi,2:3) ; roiprops(roi,4:5)]; %in image x1 y1 x2 y2   
         expinfo.driftxy=roiprops(roi,7:8); %deltax, deltay %between startframe and endframe
         expinfo.kymowidth=roiprops(roi,6); %adjust with neighbours nearby
         expinfo.channelshift=roiprops(roi,9); %manually per ROI....          
         if length(roiprops(1,:))>9, expinfo.pad_it=roiprops(roi,10); end
         
         
function expinfo=get_roiproperties_of_Atto_condensin_42kb_nicking(roi); 
        %'Atto_condensin_22kb_non_nicking\'
        expinfo.labelname='Condensin'; %no condensin channel
        expinfo.tres_pk_DNA=0.4;      %used for final peak selection; fraction of first 90% %change this for more conservative peak detection
        expinfo.tres_pk_Cnd=4;       %used for final peak selection; sigmas over noise
        roiprops=[...  %no  x1  y1  x2   y2     wd  drx dry chanshift
                        1   9  9   26   105    15  -4  -3   0; %shift upwards (-) ;leftwards (-)
                        2   19  2   28   84    14  -7  -4   0;
                        3   22  8   18  93   18   -6  -3  0;
                         4   24  8   23   153    15   -7   -2  0;
%                         5   8  6   11   99   18  -4   -3  0;
%                         6   18  4  19   88   10   -5   -2  0;
%                         7   25  23   21   102  11  -6   -3  0;
%                         8   36  13   31  93    12   -7  -2 0;
%                          9   24  10   21   59    11   -4   -3 0;
%                          10  20  12    18   48    10   -2   -1  0; %might be surface-affected
%                          11  15  11   16   140   20   0   0  0;
%                          11  15  11   16   140   20   0   0  0;
];
         expinfo.endpoints_xy=[roiprops(roi,2:3) ; roiprops(roi,4:5)]; %in image x1 y1 x2 y2   
         expinfo.driftxy=roiprops(roi,7:8); %deltax, deltay %between startframe and endframe
         expinfo.kymowidth=roiprops(roi,6); %adjust with neighbours nearby
         expinfo.channelshift=roiprops(roi,9); %manually per ROI.... 

function expinfo=get_roiproperties_of_2020_01_13_MukBEF_msd_wATP(roi); 
        %'2019_07_26 condensin_supercoil_no_ATP\'
        expinfo.labelname='MukBEF';
        expinfo.tres_pk_DNA=0.15;      %used for final peak selection; fraction of first 90%
        expinfo.tres_pk_Cnd=4;       %used for final peak selection; sigmas over noise
        roiprops=[...  %no  x1  y1  x2   y2     wd  drx dry chanshift
                        1   15  6   14   86     25   0   0  0;
                        2   10  14   40   77    77   0   0  0;
                        3   27  11   26   64    20   0   0  0;
                        4   13  12   12   87    15   0   0  0;
                        5   16  12   18   104   20   0   0  0;
                        6   11  5    18   112   20   0   0  0;
                        7   13  10   31   97    20   0   0  0;
                        8   11  10   15   77    20   0   0  0;
                        9   12  6    12   63    20   0   0  0;
                        10  15  7    10   72    18   0   0  0;
                        11  15  11   16   140   20   0   0  0;
                        ];
         expinfo.endpoints_xy=[roiprops(roi,2:3) ; roiprops(roi,4:5)]; %in image x1 y1 x2 y2   
         expinfo.driftxy=roiprops(roi,7:8); %deltax, deltay %between startframe and endframe
         expinfo.kymowidth=roiprops(roi,6); %adjust with neighbours nearby
         expinfo.channelshift=roiprops(roi,9);; %manually per ROI....  

function expinfo=get_roiproperties_of_2020_05_05_data_sc_cnd(roi); 
        %'2019_07_26 condensin_supercoil_no_ATP\'
        expinfo.labelname='Condensin';
        expinfo.tres_pk_DNA=0.15;      %used for final peak selection; fraction of first 90%
        expinfo.tres_pk_Cnd=4;       %used for final peak selection; sigmas over noise
        roiprops=[...  %no  x1  y1  x2   y2     wd  drx dry chanshift
                                              ];
         expinfo.endpoints_xy=[roiprops(roi,2:3) ; roiprops(roi,4:5)]; %in image x1 y1 x2 y2   
         expinfo.driftxy=roiprops(roi,7:8); %deltax, deltay %between startframe and endframe
         expinfo.kymowidth=roiprops(roi,6); %adjust with neighbours nearby
         expinfo.channelshift=roiprops(roi,9);; %manually per ROI....  


function expinfo=get_roiproperties_of_2019_09_02_NumberOfCondensinPerDNA(roi);
        expinfo.labelname='Condensin';
        expinfo.tres_pk_DNA=0.15;      %used for final peak selection; fraction of first 90%
        expinfo.tres_pk_Cnd=4;       %used for final peak selection; sigmas over noise
        
        switch roi
            case 1
                      expinfo.endpoints_xy=[11 26 ; 18 138]; %in image x1 y1 x2 y2   
                      expinfo.driftxy=[0 0]; %deltax, deltay %between startframe and endframe
                      expinfo.kymowidth=20; %adjust with neighbours nearby
                      expinfo.channelshift=0; %manually per ROI.... 
           case 2
                      expinfo.endpoints_xy=[23 6 ; 12 101]; %in image x1 y1 x2 y2   
                      expinfo.driftxy=[0 0]; %deltax, deltay %between startframe and endframe
                      expinfo.kymowidth=18; %adjust with neighbours nearby
                      expinfo.channelshift=0; %manually per ROI.... 
            case 3
                      expinfo.endpoints_xy=[25 5 ; 13 100]; %in image x1 y1 x2 y2   
                      expinfo.driftxy=[0 0]; %deltax, deltay %between startframe and endframe
                      expinfo.kymowidth=16; %adjust with neighbours nearby
                      expinfo.channelshift=0; %manually per ROI.... 
            case 4
                      expinfo.endpoints_xy=[28 5 ; 18 69]; %in image x1 y1 x2 y2   
                      expinfo.driftxy=[0 0]; %deltax, deltay %between startframe and endframe
                      expinfo.kymowidth=26; %adjust with neighbours nearby
                      expinfo.channelshift=0; %manually per ROI.... 
            case 5
                      expinfo.endpoints_xy=[20 2 ; 12 111]; %in image x1 y1 x2 y2   
                      expinfo.driftxy=[0 0]; %deltax, deltay %between startframe and endframe
                      expinfo.kymowidth=24; %adjust with neighbours nearby
                      expinfo.channelshift=0; %manually per ROI.... 
            case 6
                      expinfo.endpoints_xy=[22 5 ; 11 81]; %in image x1 y1 x2 y2   
                      expinfo.driftxy=[0 0]; %deltax, deltay %between startframe and endframe
                      expinfo.kymowidth=24; %adjust with neighbours nearby
                      expinfo.channelshift=0; %manually per ROI.... 
            case 7
                      expinfo.endpoints_xy=[10 5 ; 12 91]; %in image x1 y1 x2 y2   
                      expinfo.driftxy=[0 0]; %deltax, deltay %between startframe and endframe
                      expinfo.kymowidth=20; %adjust with neighbours nearby
                      expinfo.channelshift=0; %manually per ROI.... 
            case 8
                      expinfo.endpoints_xy=[18 1 ; 10 106]; %in image x1 y1 x2 y2   
                      expinfo.driftxy=[0 0]; %deltax, deltay %between startframe and endframe
                      expinfo.kymowidth=22; %adjust with neighbours nearby
                      expinfo.channelshift=0; %manually per ROI.... 
            case 9
                      expinfo.endpoints_xy=[19 9 ; 10 110]; %in image x1 y1 x2 y2   
                      expinfo.driftxy=[0 0]; %deltax, deltay %between startframe and endframe
                      expinfo.kymowidth=22; %adjust with neighbours nearby
                      expinfo.channelshift=0; %manually per ROI.... 
            case 10
                      expinfo.endpoints_xy=[24 2 ; 8 107]; %in image x1 y1 x2 y2   
                      expinfo.driftxy=[0 0]; %deltax, deltay %between startframe and endframe
                      expinfo.kymowidth=16; %adjust with neighbours nearby
                      expinfo.channelshift=0; %manually per ROI.... 
        end
        

function expinfo=get_roiproperties_of_2019_07_26_condensin_supercoil_no_ATP(roi); 
        %'2019_07_26 condensin_supercoil_no_ATP\'
        expinfo.labelname='Condensin';
        expinfo.tres_pk_DNA=0.15;      %used for final peak selection; fraction of first 90%
        expinfo.tres_pk_Cnd=4;       %used for final peak selection; sigmas over noise
        switch roi          
           case 1 
                      expinfo.endpoints_xy=[27 2 ; 13 86]; %in image x1 y1 x2 y2   
                      expinfo.driftxy=[0 0]; %deltax, deltay %between startframe and endframe
                      expinfo.kymowidth=11; %adjust with neighbours nearby
                      expinfo.channelshift=0; %manually per ROI....  
            case 2
                      expinfo.endpoints_xy=[17 13 ; 25 116]; %in image x1 y1 x2 y2   
                      expinfo.driftxy=[0 0]; %deltax, deltay %between startframe and endframe
                      expinfo.kymowidth=8; %adjust with neighbours nearby
                      expinfo.channelshift=0; %manually per ROI....  
            case 3
                      expinfo.endpoints_xy=[16 11 ; 23 118]; %in image x1 y1 x2 y2   
                      expinfo.driftxy=[0 0]; %deltax, deltay %between startframe and endframe
                      expinfo.kymowidth=11; %adjust with neighbours nearby
                      expinfo.channelshift=0; %manually per ROI....  
            case 4
                      expinfo.endpoints_xy=[15 8 ; 26 107]; %in image x1 y1 x2 y2   
                      expinfo.driftxy=[0 0]; %deltax, deltay %between startframe and endframe
                      expinfo.kymowidth=11; %adjust with neighbours nearby
                      expinfo.channelshift=0; %manually per ROI....  
            case 5
                      expinfo.endpoints_xy=[9 8 ; 22 101]; %in image x1 y1 x2 y2   
                      expinfo.driftxy=[0 0]; %deltax, deltay %between startframe and endframe
                      expinfo.kymowidth=11; %adjust with neighbours nearby
                      expinfo.channelshift=0; %manually per ROI.... 

        end 


function expinfo=get_roiproperties_of_2019_07_15_condensin_supercoil(roi);
        %'2019_07_15 condensin_supercoil\' 
        expinfo.labelname='Condensin';
        expinfo.tres_pk_DNA=0.15;      %used for final peak selection; fraction of first 90%
        expinfo.tres_pk_Cnd=4;       %used for final peak selection; sigmas over noise
        switch roi
              case 1 
                  expinfo.endpoints_xy=[18 10 ; 15 65]; %in image x1 y1 x2 y2   
                  expinfo.driftxy=[-2 0]; %deltax, deltay %between startframe and endframe
                  expinfo.kymowidth=8; %adjust with neighbours nearby
                  expinfo.channelshift=0.6022; %manually per ROI....  
              case 2 
                  expinfo.endpoints_xy=[12 7 ; 10 63]; 
                  expinfo.driftxy=[-1 0];
                  expinfo.kymowidth=8;
                  expinfo.channelshift=-1.5484; %manually per ROI....  

              case 3 
                  expinfo.endpoints_xy=[32 11; 21 101];
                  expinfo.driftxy=[-9 0];
                  expinfo.kymowidth=6;
                  expinfo.channelshift=1.4132; %manually per ROI....  

              case 4 
                  expinfo.endpoints_xy=[13 9; 14 87];
                  expinfo.driftxy=[-7 0];
                  expinfo.kymowidth=8;
                  expinfo.channelshift=1.6989; 

              case 5 
                  expinfo.endpoints_xy=[15 12; 15 95];
                  expinfo.driftxy=[-4 0];
                  expinfo.kymowidth=6;
                  expinfo.channelshift=2.0399; 
              case 6 
                  expinfo.endpoints_xy=[16 11; 18 105];
                  expinfo.driftxy=[-9 0];
                  expinfo.kymowidth=5;
                  expinfo.channelshift=2.3349; 
              case 7 
                  expinfo.endpoints_xy=[13 5; 10 98];
                  expinfo.driftxy=[-10 0];
                  expinfo.kymowidth=8;
                  expinfo.channelshift=2.5991; 
              case 8 
                  expinfo.endpoints_xy=[14 8; 16 88];
                  expinfo.driftxy=[-10 0];
                  expinfo.kymowidth=9;
                  expinfo.channelshift=1.7419; 
              case 9 
                  expinfo.endpoints_xy=[15 7; 16 99];
                  expinfo.driftxy=[-10 0];
                  expinfo.kymowidth=7;
                  expinfo.channelshift=1.7143; 
              case 10 
                  expinfo.endpoints_xy=[16 7; 16 108];
                  expinfo.driftxy=[-10 0];
                  expinfo.kymowidth=10;
                  expinfo.channelshift=1.8618; 
              case 11 
                  expinfo.endpoints_xy=[15 5; 15 87];
                  expinfo.driftxy=[-4 0];
                  expinfo.kymowidth=7;
                  expinfo.channelshift=1.2750; 
              case 12 
                  expinfo.endpoints_xy=[15 5; 16 110];
                  expinfo.driftxy=[-5 0];
                  expinfo.kymowidth=8;
                  expinfo.channelshift=0.9677; 
        end
            
