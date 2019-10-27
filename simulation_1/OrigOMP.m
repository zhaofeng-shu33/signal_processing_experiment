function [OmpResult] = OrigOMP(y,A,S)
% #########################################################################
% Orthogonal Matching Pursuit Algorithm for Sparse Signal Recovery --- OMP
% [OmpResult] = OrigOMP(y,A,S) solves the underdetermined system y=Ax where
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
% Author: MaZeqiang
% Electronic Engineering Department,Tsinghua University
% Date: 2012-12-7
% #########################################################################

%disp('Original OMP function being called') ;
%fprintf('----------------------------------\n') ;

if nargin ~= 3
    disp('error: Incorrect number of input arguments') ;
    return ;
end

[N M] = size(y) ;
[K Ns] = size(A) ;
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


LAMBDA = [] ;
r = y ;

for k = 1:S
    [~,lambda] = max(r'*A(:,1:Ns)) ;
    LAMBDA = [LAMBDA lambda] ;
    Alm = A(:,LAMBDA) ;  % the support set
    P = Alm*inv(Alm'*Alm)*Alm' ;  % projection matrix
    a = P*y ;   % orthogonal projection
    r = y -a ;  %residual error
end

Alm = A(:,LAMBDA) ;
tmp  = pinv(Alm)*a ;%Pseudoinverse.

OmpResult = zeros(Ns,1) ;
OmpResult(LAMBDA) = tmp ;

%fprintf('OMP Status: OMP solved\n') ;
%fprintf('----------------------------------\n') ;

end