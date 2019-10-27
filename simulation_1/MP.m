function [MPResult] = MP(y,A,S)
% #########################################################################
% Matching Pursuit Algorithm for Sparse Signal Recovery --- MP
% [MPResult] = MP(y,A,S) solves the underdetermined system y=Ax where
% A is a matrix of dimensions [N by Ns], y is the observation vector [N by 
% 1], x is the estimation sparse signal of [Ns by 1].  Ns>>N.
% 
% This function is only for single snapshot case, meaning the ovservation y
% is just a column vector, not a matrix.
% 
% Notation: >N: number of sensors
%           >Ns: number of sources
%           >S: sparsity degree,meaning the number of nonzeros elements in
%           the sparse vector x
%           >A: the measurement matrix
% 
% Author: ZhaoFeng
% Electronic Engineering Department,Tsinghua University
% Date: 2017-11-19
% Acckowledgement: Teaching Assistant
% #########################################################################


if nargin ~= 3
    disp('error: Incorrect number of input arguments') ;
    return ;
end

[N,M] = size(y) ;
[K,Ns] = size(A) ;
if M ~= 1
    disp('error: input arguments dimensions must agree') ;
    return ;
end

if N~= K
    disp('error: input arguments dimensions must agree') ;
    return ;
end
    
if N>=Ns
    disp('warning: A may not be an underdetermined matrix') ;
end

if S>=Ns
    disp('error: the sparse degree must not greater than colnum of A') ;
    return ;
end


LAMBDA = zeros([S,1]) ;
r = y ;
MPResult = zeros(Ns,1);
for k = 1:S
    [~,lambda] = max(r'*A(:,1:Ns)) ;
    LAMBDA(k)=lambda ;
    Alm = A(:,LAMBDA(k)) ;  % the support set
    tilde_theta=(Alm'*Alm)\(Alm'*y);
    MPResult(LAMBDA(k))=tilde_theta ;    
    r = r- Alm*tilde_theta; %residual error
end





end