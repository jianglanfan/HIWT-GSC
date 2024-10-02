% =========================================================================
% Project: HIWT-GSC
% Repository: https://github.com/jianglanfan/HIWT-GSC
%
% File Name: myLogistic.m
%
% Description:
%   This function computes the logistic regression loss and its gradient 
%   given the input parameters. It returns both the loss and the gradient 
%   if requested.
%
% Usage:
%   [f, df] = myLogistic(x, A, y)
%
% Inputs:
% - x    ---- Coefficient vector
% - A    ---- Design matrix
% - y    ---- Observed response vector
%
% Outputs:
% - f    ---- The computed logistic regression loss
% - df   ---- The gradient of the loss with respect to x (optional)
%
% References:
%   No specific references for this function.
% =========================================================================

function [f df] = myLogistic(x,A,y)
% The function computes the logistic regression loss (f) and gradient (df)
    ns = size(A,2);
    u = A'*x;
    f = sum(log(1+exp(-abs(u))) + double(u>0).*u - y.*u)/ns;
    if nargout >= 2
        df = A*(1./(1+exp(-u)) - y)/ns;
    end
