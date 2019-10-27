function [system_output,current_error_square,omega_out] = lms_system( system_input_index,get_parameter)
%LMS_SYSTEM 此处显示有关此函数的摘要
%   此处显示详细说明
global omega;
global step_size;
global system_input_sequence;
global random_sequence;
global sequence_length;
% omega is coeffients of FIR, 1 times 11 vector
% last_system_input is last 11 inputs
if nargin==2
% state is not empty
% reset
    if get_parameter==1
       omega_out=omega;
    else
        omega=zeros([1,11]);        
    end
    system_output=0;
    current_error_square=0;
else
    %output
    system_output=dot(system_input_sequence(system_input_index-10:...
        system_input_index),omega);
    current_error=random_sequence(system_input_index-7)-system_output;
    omega=omega+step_size*system_input_sequence(system_input_index-10:...
        system_input_index)*current_error;
    current_error_square=current_error*current_error;
    %if(sequence_length==system_input_index)
    %    disp(omega)
    %end
end

