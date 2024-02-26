function needtoflip=check_local_symmetry(chro_map,rel_range); 
     order=2;
     prf=max(chro_map, [], 'omitnan');
     LP=length(prf);
     mid=round(LP/2);
     mask_prf=(OnePeak(1:LP,mid,rel_range*LP,'Periodic')).^order;     
     masked=prf.*mask_prf;    
     leftsum=sum(masked(1:mid-1));
     rightsum=sum(masked(mid:end));
     localsymmetry=(rightsum-leftsum)/sum(masked);
     needtoflip=(localsymmetry<0);
     if 0
         subplot(2,1,1);
            pcolor(chro_map); shading flat; colormap jet;
            
         subplot(2,1,2);
            plot(prf); hold on;
            plot(masked, 'r-'); hold on;   
            plot([mid mid], [0 max(masked)], 'r--');
            legend('profile', 'masked profile', 'parB focus');
            xlabel('position');
            ylabel('signal');
            [~]=ginput(1);
     end
     