 function [ph_pic, msk_pic, c2_pic, c3_pic, c4_pic]=get_channel_ID(channel_stack, initval)
    %allocate pictures for further analysis
    %pictures were loaded by name to stack as [ma c1 c2 c3 c4]
    %default chan_ID=[1 0 2 3 1];  %ph ms ch or tr 
    ID=initval.chan_ID; 
    ph_pic=channel_stack(:,:,ID(1)+1);
    msk_pic=channel_stack(:,:,ID(2)+1);
    c2_pic=channel_stack(:,:,ID(3)+1); %ch
    c3_pic=channel_stack(:,:,ID(4)+1); %ori
    c4_pic=channel_stack(:,:,ID(5)+1); %ter 