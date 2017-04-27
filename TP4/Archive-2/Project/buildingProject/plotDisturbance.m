function [] = plotDisturbance( refDist )
    timeLine = linspace(0, 8, 576);
    figure();
    subplot(3,1,1);
    plot(timeLine,refDist(1,:)');
    xlabel('Days');
    ylabel('Temperature (°C)');
    title('Disturbance input');

    subplot(3,1,2);
    plot(timeLine,refDist(2,:)');
    xlabel('Days');
    ylabel('Solar gains (kW)');
    
    subplot(3,1,3);
    plot(timeLine,refDist(3,:)');
    xlabel('Days');
    ylabel('Internal gains (kW)');


end

