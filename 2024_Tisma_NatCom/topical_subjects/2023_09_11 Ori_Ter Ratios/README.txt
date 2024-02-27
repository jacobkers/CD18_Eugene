## Ori-ter ratios project
We started with Oufti data. Issue here was the poor quality of the ter channel.
Therefore, I tried another approach, 'pek-peeling' where we re-analyze the images and try to get a reasoanble treshold for spots to accept. This did not give a major improvement.
Thus, we went back to the oufti data and in the end did a bit more critical pre-screening of the cells: we would only accept spots that contain a sufficient fraction of the total signal.
This step would reject cells with a very homogenous fluorescence, which would give an unreasonable number of ter spots.
Judging by eye, and holding the same approach for the ori channel (which was cleaner) we chose 0.15 as content fraction limit.





## Info Tisma:----------------------------------

So I have made the cell outlines in Oufti and the spot detection for 0-30-60-120-180 min replication halt.

The 180min still needs to be refined but I won't have time so as long as the code runs it can be repeated later too.

The spot locations and numbers for ter are found in cellList.meshData{n1,n2}{m1,m2}.spots.
the corresponding for ori are found in cellList.meshData{n1,n2}{m1,m2}.spots2

.mat files are found in M:\tnw\bn\cd\Shared\Tisma\Bacillus subtilis chromosome imaging\Widefield_and_deconvolution\20230906_BSG5522_IPTG_timecourse
always use the file that ends up with _02 (e.g., BSG5522_60min_c1_Stack_SpotDetection_02.mat)

So in here we basically needs:
1) number of ori per cell (can be like plots at the top of Figure 2)
2) number of ter per cell
3) average ori/ter ratio in each condition 

I will be in a conference almost entire next week but I am available by mail for any discussion about this! I will also try to finalize the 180min condition but for now that might also be done with this preliminary stuff.

Best,

Tisma

Miloš Tišma
Office : F0.170
Delft University of Technology
Kavli Institute of Nanoscience Delft
Department of Bionanoscience
Van der Maasweg 9
2629 HZ Delft
The Netherlands

