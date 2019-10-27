function OMPExperiment1()
global k0
[mse,hat_x,x]=MPOMP_wrapper(4*k0,10);
fprintf('experiment_1 error: %E,\n',mse)   
%plot...
stem(x)
hold on
stem(hat_x,'Marker','*','MarkerSize',6)
legend('Original Signal','CS with OMP')
xlim([1,1024])
ylim([0,250])
title('Compress Sensing')
saveas(gcf,'Compress_Sensing_1','epsc')
end