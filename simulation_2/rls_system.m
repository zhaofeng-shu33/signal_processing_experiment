function [system_output,current_error_square] = rls_system( system_input_index,~ )
global omega;
global P;
global system_input_sequence;
global random_sequence;
global forget_parameter;
% omega is coeffients of FIR, 1 times 11 vector
% last_system_input is last 11 inputs
if nargin==2
% state is not empty
% reset
    omega=zeros([1,11]);
    P=10*eye(11);
else
    %output
    xn=system_input_sequence(system_input_index-10:system_input_index);
    dn=random_sequence(system_input_index-7);
    system_output=dot(xn,omega);    
    current_error=dn-system_output;    
    current_error_square=current_error*current_error;    
    kn=P*xn'/(forget_parameter+xn*P*xn');
    epsilon_n=dn-omega*xn';
    omega=omega+epsilon_n*kn';
end

