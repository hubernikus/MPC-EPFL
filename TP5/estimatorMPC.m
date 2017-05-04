function [x_hat2, d_hat2] estimatorMPC(A_err, B_err, L, C_err, u, y, x_hat, d_hat)
%ESTIMATORMPC - estimation of the funciton

dh_hat2 = A_err*[x_hat; d_hat]+B_err*u+L*[C_err,y];
x_hat2 = dh_hat2(1:2);
d_hat2 = dh_hat2(3);