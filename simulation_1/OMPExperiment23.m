function OMPExperiment23(stage)
global k0
M_seq=floor(linspace(k0,5*k0,16));
error_list=zeros([1,length(M_seq)]);
for j=1:length(M_seq)
  if(stage==2)
  [mse,~,~]=MPOMP_wrapper(M_seq(j));
  else
  [mse,~,~]=MPOMP_wrapper(M_seq(j),20,'MP');
  end
  error_list(j)=mse;
  disp(j);
end
%plot...
figure(stage)
semilogy(M_seq,error_list,'b')
hold on
semilogy(M_seq,error_list,'b*')
xlabel('M')
ylabel('Mean Square Error')
title('Compressed Sensing')
if(stage==2)
saveas(gcf,'Compress_Sensing_2','epsc')
else
saveas(gcf,'Compress_Sensing_3','epsc')
end    
end

