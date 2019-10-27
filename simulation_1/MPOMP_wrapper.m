function [mse,hat_x,x] = MPOMP_wrapper(M,repeat_time,method)
%MPOMP_WRAPPER 此处显示有关此函数的摘要
%   此处显示详细说明
global N k0 random_data_x
OMP_method=true;
if nargin == 1
    repeat_time=20;
    disp('method: OMP, repeat time:20\n') ;
elseif nargin == 2
    assert(repeat_time>=1)
    fprintf('method: OMP, repeat time:%d\n',repeat_time);
else    
    assert(repeat_time>=1)
    OMP_method=strcmp(method,'MP');
    assert(strcmp(method,'OMP')==0 || strcmp(method,'MP')==0)
    fprintf('method: %s, repeat time:%d\n',method,repeat_time);    
end
%initialization
error=0;
for i=1:repeat_time
Phi=randn([M,N])/sqrt(M);
x=zeros([N,1]);
x(randsample(N,k0))=random_data_x;
y=Phi*x;
if(OMP_method)
hat_x=OMP(y,Phi,k0);
else %MP
hat_x=MP(y,Phi,k0);
end
error_increment=norm(x-hat_x)^2/norm(x)^2;
%disp(error_increment)
error=error+error_increment;
end
mse=error/repeat_time;
end