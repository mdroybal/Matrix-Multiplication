%Monte Roybal
%CS_577 Parallel and Distributed Programming
%5-10-2018
%Dr. Gil Gallegos
%Freemat N x M Final

clc
clear

A = [4 -1 0 0;
    -1 4 -1 0;
    0 -1 4 -1;
    0 0 -1 4;
    0 -2 0 -2;];
    
b = [10 0 50 20 30]';
    
A_T = A';

A_star = A_T * A;

b_star = A_T * b;

%x = (A_star)_-1 * b_star

A_star_inv = inv(A_star);

x = A_star_inv * b_star

%Moore-Penrose Pseudo-Inverse

A_mp_inv = pinv(A_star);

x_mp = A_mp_inv * b_star

%A transpose using for loop

size_A = size(A);
N = size_A(1);
M = size_A(2);

for i=1:M
    for j=1:N
        A_T(i,j) = A(j,i);
    end
end

B = A_T;

%Matrix Matrix Multiply
sum_MM = 0;
for i=1:M
    for j=1:M
        for k=1:N
            sum_MM = sum_MM + B(i,k) * A(k,j);
        end
        C(i,j) = sum_MM;
        sum_MM = 0;
    end
end

C







