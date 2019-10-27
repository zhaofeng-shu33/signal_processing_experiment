clc
clear all
close all


m=3;%目标个数
f_ms=[-28 5 25].*pi/180 ; %0.3*(1:m)-0.66;%真实接收角
% N_s=100;%划分空间网格数
SNR=20;%信噪比
N=10;%采样数
Q=2;
snap=1;
tic

f_m_s= [-50:2:50]*pi/180  ; %-0.5:(1/(N_s-1)):0.5;%划分接收角
S_m_s=exp(j*((1:N)-1)'*2*pi*f_m_s);%接收角字典
Fai=S_m_s;%总体字典(内存可能不足！！！！！！！！！！！！)
toc


S_ms=exp(j*((1:N)-1)'*2*pi*f_ms);%目标接收角导引矢量

  %产生目标回波
Y=10^(SNR/20)*S_ms*randn(m,snap)+randn(N,snap);



%% 原始重构算法
tic
Labda=OrigOMP(Y,Fai,m);
L_1=find(Labda~=0);%重构扩充位置

toc


figure(3)
plot(f_ms*180/pi,'r.','markersize',30)
hold on
plot(f_m_s(L_1)*180/pi,'o','markersize',10)
grid on
legend('true targets','estimated targets',1)
xlabel('目标')
ylabel('目标角度(度)')
% axis([0 4 -0.5 0.5])


