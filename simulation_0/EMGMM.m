function EMGMM()
%parameter initialization
error_report_counter=1;
mu=[10,3;1,1;5,4];
sigma=[1,1;1.5,1.5;2,2];
C=cat(3,[1,0;0,1],[1.5,0;0,1.5],[2,0;0,2]);

%generate 100 N(mu_1,C_1) 
%generate 100 N(mu_2,C_2) 
%generate 100 N(mu_3,C_3)
%sample container: Y=[], len(Y)=300, Y(i)=[x,y]
% C_i is decoupled, generate each component respectively

%random number generation
Y=zeros([300,2]);
for i=1:3
    for j=1:2
        for k=1:100
            Y(100*(i-1)+k,j)=random('Normal',mu(i,j),sigma(i,j));
        end
    end
end

%GMM plotting routine
if(nargin==1)
scatter(Y(:,1),Y(:,2),'X')
title('sample distribution')
saveas(gcf,'sample_distribution','epsc')
return;
end

%generate 3 datasets
w_initial=[0.3,0.3,0.4];
mu_initial=mu-random('normal',0,1,3,2);
C_initial=C;
y1=EMGMM_task(w_initial,mu_initial,C_initial);
w_initial=[0.25,0.25,0.5];
mu_initial=mu-random('normal',0,1.5,3,2);
C_initial=C;
y2=EMGMM_task(w_initial,mu_initial,C_initial);
w_initial=[0.45,0.4,0.15];
mu_initial=mu-random('normal',0,2.5,3,2);
C_initial=cat(3,[1,0;0,1],[1,0;0,1],[1,0;0,1]);
y3=EMGMM_task(w_initial,mu_initial,C_initial);
clf
plot(1:length(y1),y1,'r','LineWidth',2);
hold on
plot(1:length(y2),y2,'b','LineWidth',2)
plot(1:length(y3),y3,'k','LineWidth',2)
legend('small','medium','big')
scatter(1:length(y1),y1,'r')
scatter(1:length(y2),y2,'b')
scatter(1:length(y3),y3,'k')
title('influence of initial values on EM algorithm')
xlabel('iteration time')
ylabel('likelyhood function value')
saveas(gcf,'likelyhood_funtion_change','epsc')
function y=EMGMM_task(w_initial,mu_initial,C_initial)%return LikelyhoodFunctionAverage
mMax=500;
LikelyhoodFunction=Inf([mMax,1]);
%used for estimating the final error
wLast=w_initial;
muLast=mu_initial;
CLast=C_initial;

epsilon_estimation=0.001;
%average over different random samples
EMGMM_inner(w_initial,mu_initial,C_initial,epsilon_estimation);

end_pointer=1;
while(LikelyhoodFunction(end_pointer)~=Inf)    
    end_pointer=end_pointer+1;
end
end_pointer=end_pointer-1;
average_times=1;
LikelyhoodFunctionAverage=zeros([average_times,end_pointer]);
for it=1:average_times
    LikelyhoodFunctionAverage(it,:)=LikelyhoodFunction(1:end_pointer)';
    if(it~=average_times)
        EMGMM_inner(w_initial,mu_initial,C_initial,epsilon_estimation);    
    end
end
%disp the initial error and the estimated error
error_initial=sprintf('%d & %0.2f & %0.2f & %0.2f & %0.2f & %0.2f & %0.2f \\\\',error_report_counter,norm(w_initial-[1,1,1]/3),...+
norm(wLast-[1,1,1]/3),norm(mu_initial-mu),norm(muLast-mu),...+
norm(C_initial(:,:,1)-C(:,:,1))+norm(C_initial(:,:,2)-C(:,:,2))+norm(C_initial(:,:,3)-C(:,:,3)),...+
norm(CLast(:,:,1)-C(:,:,1))+norm(CLast(:,:,2)-C(:,:,2))+norm(CLast(:,:,3)-C(:,:,3)));
disp(error_initial);
error_report_counter=error_report_counter+1;
%then return LikelyhoodFunctionAverage
y=mean(LikelyhoodFunctionAverage,1);

function EMGMM_inner(given_w,given_mu,given_C,given_epsilon)
%EM Algorithm Implementation
%S1
%set maximum iteration times
K=3;
N=300;
%m: iteration pointer
m=1;
epsilon=given_epsilon;%0.01;
gamma=zeros([N,K]);
smallN=zeros([1,K]);
%we only save the current w,mu,C and last w,mu,C
%initial guess of w,mu,C
wLast=given_w;%[0.33,0.33,0.34];
muLast=given_mu;
CLast=given_C;
wCurrent=wLast;
muCurrent=muLast;
CCurrent=CLast;
LikelyhoodFunction(1)=CalculateLikelyhood(wLast,muLast,CLast);
while(m<mMax)
%S2:E
    for i=1:N
        tmp_sum=0;
        for k=1:K
            gamma(i,k)=wLast(k)*CalculateNormalPDF(Y(i,:),muLast(k,:),CLast(:,:,k));
            tmp_sum=tmp_sum+gamma(i,k);
        end
        for k=1:K
            gamma(i,k)=gamma(i,k)/tmp_sum;
        end
    end
    for k=1:K
        smallN(k)=sum(gamma(:,k));
    end

%S3:M
    wCurrent=smallN/N;
    for k=1:K    
    muCurrent(k,:)=[0,0];
        for i=1:N
            muCurrent(k,:)=muCurrent(k,:)+gamma(i,k)*Y(i,:);
        end
    muCurrent(k,:)=muCurrent(k,:)/smallN(k);
    end
    for k=1:K    
    CCurrent(:,:,k)=[0,0;0,0];
        for i=1:N
            CCurrent(:,:,k)=CCurrent(:,:,k)+gamma(i,k)*(Y(i,:)-muCurrent(k,:))'*(Y(i,:)-muCurrent(k,:));
        end
    CCurrent(:,:,k)=CCurrent(:,:,k)/smallN(k);
    end

%S4
    LikelyhoodFunction(m+1)=CalculateLikelyhood(wCurrent,muCurrent,CCurrent);
    assert(LikelyhoodFunction(m+1)>LikelyhoodFunction(m));
    if(LikelyhoodFunction(m+1)<LikelyhoodFunction(m)+epsilon)
        break
    else
        m=m+1;
        wLast=wCurrent;
        muLast=muCurrent;
        CLast=CCurrent;
    end
end
function y=CalculateNormalPDF(y_ii,mu_ii,C_ii)
%dimension=2
%y_ii,mu_ii are 2d row vectors
%C_ii is 2\times 2 matrix
    tmp_diff=y_ii-mu_ii;
    y=(1/(2*pi*sqrt(det(C_ii))))*exp(-tmp_diff*(C_ii\tmp_diff')/2);
end
function L=CalculateLikelyhood(w_i,mu_i,C_i)
    L=0;
    for iii=1:N
        L_inner=0;
        for jjj=1:K
            L_inner=L_inner+w_i(jjj)*CalculateNormalPDF(Y(iii,:),mu_i(jjj,:),C_i(:,:,jjj));
        end
        L=L+log(L_inner);
    end
    L=L/N;
end
end
end
end