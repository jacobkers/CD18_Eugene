function expinfo=A002_EK_Condensin_with_plectonemes_expinfo(expi,roi);
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
    case 1,   expinfo=get_roiproperties_of_2019_07_15_condensin_supercoil(roi);          
    case 2,   expinfo=get_roiproperties_of_2019_07_26_condensin_supercoil_no_ATP(roi);         
    case 3,   expinfo=get_roiproperties_of_2019_09_02_NumberOfCondensinPerDNA(roi);  
    case 4,   expinfo=get_roiproperties_of_2020_01_13_MukBEF_msd_wATP(roi);    
    case 5,   expinfo=get_roiproperties_of_2019_7_26_condensin_supercoil_with_ATP(roi);  
    case 6,   expinfo=get_roiproperties_of_WT_condensin_22kb_non_nicking(roi); 
    case 7,   expinfo=get_roiproperties_of_Atto_condensin_22kb_non_nicking(roi); 
    case 8,   expinfo=get_roiproperties_of_42kb_sc_control(roi); 
end
function expinfo=get_roiproperties_of_42kb_sc_control(roi); 
        %'Atto_condensin_22kb_non_nicking\'
        expinfo.labelname='DNA'; %no condensin channel
        expinfo.tres_pk_DNA=0.4;      %used for final peak selection; fraction of first 90% %change this for more conservative peak detection
        expinfo.tres_pk_Cnd=0.4;      % no condensin channel
       % expinfo.tres_pk_Cnd=3;       %used for final peak selection; sigmas over noise
        roiprops=[...  %no  x1  y1  x2   y2     wd  drx dry chanshift
                        1   28  11   30   80    25  -8  -1   0; %shift upwards (-)
                        2   27  19   21   76    21  -2  -1   0;
                        3   20  16   12   75    17   -2  -1  0;
                        4   20  21   22   75    12   -2   0  0;
                         5   22  21   21   74   18  -1   -1  0;
                          6   17  12   16   78   10   -3   -2  0;
%                          7   31  11   17   66    12   -4   -2  0;
%                          8   11  31   22   69    12   -4  -2 0;
%                          9   24  10   21   59    11   -4   -3 0;
%                          10  20  12    18   48    10   -2   -1  0; %might be surface-affected
%                          11  15  11   16   140   20   0   0  0;
%                          11  15  11   16   140   20   0   0  0;
];
         expinfo.endpoints_xy=[roiprops(roi,2:3) ; roiprops(roi,4:5)]; %in image x1 y1 x2 y2   
         expinfo.driftxy=roiprops(roi,7:8); %deltax, deltay %between startframe and endframe
         expinfo.kymowidth=roiprops(roi,6); %adjust with neighbours nearby
         expinfo.channelshift=roiprops(roi,9); %manually per ROI....  
function expinfo=get_roiproperties_of_Atto_condensin_22kb_non_nicking(roi); 
        %'Atto_condensin_22kb_non_nicking\'
        expinfo.labelname='Condensin';
        expinfo.tres_pk_DNA=0.4;      %used for final peak selection; fraction of first 90% %change this for more conservative peak detection
        expinfo.tres_pk_Cnd=3;       %used for final peak selection; sigmas over noise
        roiprops=[...  %no  x1  y1  x2   y2     wd  drx dry chanshift
                        1   14  20   11   79    15  -2  -1   0; %shift upwards (-)
                        2   16  21   18   85    15  -3  -2   0;
                        3   17  14   17   82    17   -4  -2  0;
                        4   16  16   9   77    17   -4   -1  0;
                         5   10  15   18   80   18  -3   -2  0;
                         6   17  12   16   78   10   -3   -2  0;
                         7   31  11   17   66    12   -4   -2  0;
                         8   11  31   22   69    12   -4  -2 0;
                         9   24  10   21   59    11   -4   -3 0;
                         10  20  12    18   48    10   -2   -1  0; %might be surface-affected
%                         11  15  11   16   140   20   0   0  0;
];
         expinfo.endpoints_xy=[roiprops(roi,2:3) ; roiprops(roi,4:5)]; %in image x1 y1 x2 y2   
         expinfo.driftxy=roiprops(roi,7:8); %deltax, deltay %between startframe and endframe
         expinfo.kymowidth=roiprops(roi,6); %adjust with neighbours nearby
         expinfo.channelshift=roiprops(roi,9); %manually per ROI....  
function expinfo=get_roiproperties_of_WT_condensin_22kb_non_nicking(roi); 
        %'WT_condensin_22kb_non_nicking\'
       % expinfo.labelname='Condensin';
        expinfo.labelname=0;
        expinfo.tres_pk_DNA=0.10;      %used for final peak selection; fraction of first 90%
        expinfo.tres_pk_Cnd=4;       %used for final peak selection; sigmas over noise
        roiprops=[...  %no  x1  y1  x2   y2     wd  drx dry chanshift
                        1   25  20   28   57    20  3  0   0;
                        2   20  32   48   33    20  2  0   0;
                        3   17   8   10   45    12   0  0  0;
%                         4   13  12   14   114    18   0   0  0;
%                         5   18  15   20   105   18   0   0  0;
%                         6   11  5    18   112   15   0   0  0;
%                         7   13  10   31   97    20   0   0  0;
%                         8   11  10   15   77    20   0   0  0;
%                         9   12  6    12   63    18   0   0  0;
%                         10  15  7    10   72    18   0   0  0;
%                         11  15  11   16   140   20   0   0  0;
];
         expinfo.endpoints_xy=[roiprops(roi,2:3) ; roiprops(roi,4:5)]; %in image x1 y1 x2 y2   
         expinfo.driftxy=roiprops(roi,7:8); %deltax, deltay %between startframe and endframe
         expinfo.kymowidth=roiprops(roi,6); %adjust with neighbours nearby
         expinfo.channelshift=roiprops(roi,9); %manually per ROI....  

function expinfo=get_roiproperties_of_2020_01_13_MukBEF_msd_wATP(roi); 
        %'2019_07_26 condensin_supercoil_no_ATP\'
        expinfo.labelname='Condensin';
        expinfo.tres_pk_DNA=0.15;      %used for final peak selection; fraction of first 90%
        expinfo.tres_pk_Cnd=4;       %used for final peak selection; sigmas over noise
        roiprops=[...  %no  x1  y1  x2   y2     wd  drx dry chanshift
                        1   15  6   14   86     25   0   0  0;
                        2   13  16   40   77    15   0   0  0;
                        3   27  11   26   64    20   0   0  0;
                        4   13  12   14   114    18   0   0  0;
                        5   18  15   20   105   18   0   0  0;
                        6   11  5    18   112   15   0   0  0;
                        7   13  10   31   97    20   0   0  0;
                        8   11  10   15   77    20   0   0  0;
                        9   12  6    12   63    18   0   0  0;
                        10  15  7    10   72    18   0   0  0;
                        11  15  11   16   140   20   0   0  0;
                        ];
         expinfo.endpoints_xy=[roiprops(roi,2:3) ; roiprops(roi,4:5)]; %in image x1 y1 x2 y2   
         expinfo.driftxy=roiprops(roi,7:8); %deltax, deltay %between startframe and endframe
         expinfo.kymowidth=roiprops(roi,6); %adjust with neighbours nearby
         expinfo.channelshift=roiprops(roi,9);; %manually per ROI....  


function expinfo=get_roiproperties_of_2019_09_02_NumberOfCondensinPerDNA(roi);
        expinfo.labelname='Condensin';
        expinfo.tres_pk_DNA=0.10;      %used for final peak selection; fraction of first 90%
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
        expinfo.tres_pk_DNA=0.2;      %used for final peak selection; fraction of first 90%
        expinfo.tres_pk_Cnd=5;       %used for final peak selection; sigmas over noise
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
            
function expinfo=get_roiproperties_of_2019_7_26_condensin_supercoil_with_ATP(roi);  
        %'2019_07_15 condensin_supercoil\' 
        expinfo.labelname='Condensin';
        expinfo.tres_pk_DNA=0.15;      %used for final peak selection; fraction of first 90%
        expinfo.tres_pk_Cnd=3;       %used for final peak selection; sigmas over noise
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

%               case 4 
%                   expinfo.endpoints_xy=[13 9; 14 87];
%                   expinfo.driftxy=[-7 0];
%                   expinfo.kymowidth=8;
%                   expinfo.channelshift=1.6989; 
% 
%               case 5 
%                   expinfo.endpoints_xy=[15 12; 15 95];
%                   expinfo.driftxy=[-4 0];
%                   expinfo.kymowidth=6;
%                   expinfo.channelshift=2.0399; 
%               case 6 
%                   expinfo.endpoints_xy=[16 11; 18 105];
%                   expinfo.driftxy=[-9 0];
%                   expinfo.kymowidth=5;
%                   expinfo.channelshift=2.3349; 
%               case 7 
%                   expinfo.endpoints_xy=[13 5; 10 98];
%                   expinfo.driftxy=[-10 0];
%                   expinfo.kymowidth=8;
%                   expinfo.channelshift=2.5991; 
%               case 8 
%                   expinfo.endpoints_xy=[14 8; 16 88];
%                   expinfo.driftxy=[-10 0];
%                   expinfo.kymowidth=9;
%                   expinfo.channelshift=1.7419; 
%               case 9 
%                   expinfo.endpoints_xy=[15 7; 16 99];
%                   expinfo.driftxy=[-10 0];
%                   expinfo.kymowidth=7;
%                   expinfo.channelshift=1.7143; 
%               case 10 
%                   expinfo.endpoints_xy=[16 7; 16 108];
%                   expinfo.driftxy=[-10 0];
%                   expinfo.kymowidth=10;
%                   expinfo.channelshift=1.8618; 
%               case 11 
%                   expinfo.endpoints_xy=[15 5; 15 87];
%                   expinfo.driftxy=[-4 0];
%                   expinfo.kymowidth=7;
%                   expinfo.channelshift=1.2750; 
%               case 12 
%                   expinfo.endpoints_xy=[15 5; 16 110];
%                   expinfo.driftxy=[-5 0];
%                   expinfo.kymowidth=8;
%                   expinfo.channelshift=0.9677; 
        end
        
            
