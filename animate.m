figure(101),scatter(positions(:,1),positions(:,2),50,r(la,:),'filled')
set(gca,'FontSize',18)
    text(0.1,0.95*max(max(positions)),['t = ' num2str(t(la))],'FontSize',24)
    %text(0,0.55*max(max(positions)),['Activity = ' num2str(sum(r(la,:)))],'FontSize',24)    
    %title(['Network activity: firing rates at N=' num2str(N) ' real-space positions'])
    title(['Network activity: firing rates (N=' num2str(N) ')'])
    xlabel('Neuron x-position [mm]')
    ylabel('Neuron y-position [mm]')
    %xlabel('Neuron x-position [multiples of interneuron separation l]')
    %ylabel('Neuron y-position [multiples of interneuron separation l]')
    caxis manual
    caxis([0 4])
    %caxis([0 10]);
    %caxis([0 a]);
    %colormap jet;
    colormap(flipud(autumn))
    set(gca,'Color',[1 1 0])
    colorbar
    drawnow;
    %pause(0.5)
    %viscircles([positions(1,1) positions(1,2)],restr,'Color','k');
    xlim([0 5])
    ylim([0 5])
    set(gca,'FontSize',18)
    
    %figure(101);
    %frame = getframe(gcf);
    %writeVideo(v,frame);
