%direction of arrival
%three target
target_degree=[-28,5,25];
target=target_degree*pi/180;
signal=[1,1,1];
dLambda=1;
n=10;
noise_variance=0.1;
n_seq=1:n;
%split of search theta space
theta_M=(-50:2:50)*pi/180;
%construct H matrix
H=zeros([n,length(theta_M)]);
for i=1:length(theta_M)
    H(:,i)=exp(-1j*2*pi*(n_seq-1)*sin(theta_M(i))/dLambda);
end
%generate y through random sampling
A=zeros([n,length(target)]);
for i=1:length(target)
    A(:,i)=exp(-1j*2*pi*(n_seq-1)*sin(target(i))/dLambda);
end
observation=A*signal'+randn([n,1])*noise_variance;
%now run OMP algorithm
OmpResult=OrigOMP(observation,H,length(target));
estimated_pos=zeros([1,length(target)]);
cnt=1;
for i=1:length(OmpResult)
    if(abs(OmpResult(i))>0.1)
        estimated_pos(cnt)=-50+2*(i-1);
        cnt=cnt+1;
    end
end
assert(cnt==1+length(target));
plot(target_degree,'r.','markersize',30)
hold on
plot(estimated_pos,'o','markersize',10)
grid on
legend('true targets','estimated targets')
xlabel('目标')
ylabel('目标角度(度)')
% axis([0 4 -0.5 0.5])
