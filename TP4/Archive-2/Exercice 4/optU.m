function [flag, u] = optU(H, h, G, g, T, t)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    [zopt, fval, flag] = quadprog(H, h, G, g, T, t);
    u = zopt((N*size(A,2) + 1):(N*size(A,2) + size(B,2)));
end

