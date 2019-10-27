clc
clear all
close all


m=3;%Ŀ�����
f_ms=[-28 5 25].*pi/180 ; %0.3*(1:m)-0.66;%��ʵ���ս�
% N_s=100;%���ֿռ�������
SNR=20;%�����
N=10;%������
Q=2;
snap=1;
tic

f_m_s= [-50:2:50]*pi/180  ; %-0.5:(1/(N_s-1)):0.5;%���ֽ��ս�
S_m_s=exp(j*((1:N)-1)'*2*pi*f_m_s);%���ս��ֵ�
Fai=S_m_s;%�����ֵ�(�ڴ���ܲ��㣡����������������������)
toc


S_ms=exp(j*((1:N)-1)'*2*pi*f_ms);%Ŀ����սǵ���ʸ��

  %����Ŀ��ز�
Y=10^(SNR/20)*S_ms*randn(m,snap)+randn(N,snap);



%% ԭʼ�ع��㷨
tic
Labda=OrigOMP(Y,Fai,m);
L_1=find(Labda~=0);%�ع�����λ��

toc


figure(3)
plot(f_ms*180/pi,'r.','markersize',30)
hold on
plot(f_m_s(L_1)*180/pi,'o','markersize',10)
grid on
legend('true targets','estimated targets',1)
xlabel('Ŀ��')
ylabel('Ŀ��Ƕ�(��)')
% axis([0 4 -0.5 0.5])


