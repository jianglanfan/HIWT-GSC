% =========================================================================
% Project: HIWT-GSC
% Repository: https://github.com/jianglanfan/HIWT-GSC
%
% File Name: myLinear.m
%
% Description:
%   This function computes the linear least squares cost and its gradient 
%   given the input parameters. It returns both the cost and the gradient 
%   if requested.
%
% Usage:
%   [f, df] = myLinear(x, A, y)
%
% Inputs:
% - x    ---- Coefficient vector
% - A    ---- Design matrix
% - y    ---- Observed response vector
%
% Outputs:
% - f    ---- The computed cost (objective function value)
% - df   ---- The gradient of the cost with respect to x (optional)
%
% References:
%   No specific references for this function.
% =========================================================================

function [f df] = myLinear(x,A,y)
    u = A*x;
    f = 0.5*sum((u-y).^2);
    if nargout >= 2
        df = A'*(u-y);
    end
    