function kymo_out=kym_pre_cook_kymo_DNA(kymo);
      %clean up: background slope, bleaching
      %1 flatten background left-right
      
      prf_med=median(kymo);
      allminval=min(prf_med);
      [ff,Lp]=size(kymo);
      xax=1:Lp;
        
      [lo_L,ixL]=min(prf_med(1:ceil(Lp/2)));
      [lo_R,ixR]=min(prf_med(ceil(Lp/2):end)); ixR=ixR+ceil(Lp/2)-1;
      slopeline=polyval(polyfit([ixL ixR],[lo_L lo_R],1),xax);
      slopeplane=repmat(slopeline,ff,1);
      kymo=kymo-slopeplane+allminval;
      
      
      %2 crude bleach correct
      tax=1:ff;
      bleachline=mean([kymo(1:ff,1) kymo(1:ff,end)]');
      bleachline_nrm=bleachline/mean(bleachline);
      bleachslope=polyval(polyfit(tax,bleachline_nrm,1),tax);
      bleachplane=repmat(bleachslope',1,Lp);
      
      kymo=kymo./bleachplane;
      kymo_out=matrix_blur_jk(kymo,3);
      
      for ii=1:ff
          kymo_out(ii,:)=kymo_out(ii,:)-min(kymo_out(ii,:));
      end
      
      
      if 0          
        subplot(1,2,1);
            plot(tax(2:end),bleachline_nrm(2:end)); hold on;
            plot(tax(2:end),bleachslope(2:end),'r');
            title('bleach correction');
            legend('mean edges, norm.','linear fit');
            xlabel('time, frames');
            ylabel('norm.intensity, a.u.');
            xlim([1 ff]);
        subplot(1,2,2);           
            plot(prf_med,'b-'); hold on;
            plot(slopeline,'b-');
            plot(prf_med-slopeline+allminval,'r'); 
            title('flattening');
            legend('median','minima fit','corrected');
            xlabel('position, pixel units');
            ylabel('fluorescence intensity, a.u.');
            pause(0.5);        
            [~]=ginput(1);
            hold off
      end