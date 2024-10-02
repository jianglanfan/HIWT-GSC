% =========================================================================
% Project: HIWT-GSC
% Repository: https://github.com/jianglanfan/HIWT-GSC
%
% File Name: demo_LogisticRegression.m
% Author: Lanfan Jiang
% Date Created: September 2024
% Last Modified: September 29, 2024
% Version: N/A
%
% Description:
%   This example demonstrates how to use the HIWT-GSC for a 
%   group-sparsity-constrained logistic regression problem. 
%   We evaluate the performance of HIWT_GSC in terms of group feature
%   selection accuracy and prediction error rate. The evaluation metrics 
%   include F1-score for accuracy and ER (Error Rate) for prediction error.
%
% Usage:
%   To run the demo, simply execute the corresponding script in MATLAB:
%   >> demo_LogisticRegression
%
% License:
%   CC BY-NC 4.0
%
% References:
%   If you use this code, please cite the following paper:
%   
%   L. Jiang, Z. Huang, Y. Chen, and W. Zhu, 
%   "Iterative-Weighted Thresholding Method for 
%   Group-Sparsity-Constrained Optimization with Applications," 
%   IEEE Transactions on Neural Networks and Learning Systems, 
%   early access, 2024. 
%   DOI: 10.1109/TNNLS.2024.3454070
% =========================================================================

clear all
clc
close all
addpath(genpath(fileparts(mfilename('fullpath'))));
FIELDS={'HIWT_GSC'};        % all available solvers; 

%% Data settings 
SmallDim = true;            % whether to test with small dimensional data
if SmallDim
    p = 1000;               % signal length  
    n = 250;                % number of samples, sampling rate = 25%
    num_groups = 100;       % number of groups in the signal 
    SparseLev = 5:1:15;     % sparse level: 5%-15%
else
    p = 10000;              % signal length
    n = 5000;               % number of samples, sampling rate = 50%
    num_groups = 1000;      % number of groups in the signal 
    SparseLev = 50:10:150;  % sparse level: 5%-15%
end
gs = p/num_groups;          % number of elements in each group
nmethod = numel(FIELDS);    % number of methods involved in testing
maxnumtest = 100;           % number of tests performed for each parameter setting
len = length(SparseLev);

%% Define arrays to save the results of each test
Scputime = zeros(maxnumtest,len,nmethod);
ER = zeros(maxnumtest,len,nmethod);     % error rate 
F1 = zeros(maxnumtest,len,nmethod);     % F1-score 

dopts.seednum  = 0;         % seed number (default 0) 
rng('default'); 
noise = false;
fig = 1;
fid = 1;
printf = @(varargin) fprintf(fid,varargin{:});

for k = 1:len
    num_nz_groups = SparseLev(k);
        for l = 1:maxnumtest
            printf('\n Sampling Rate = : %g, Sparse Level = %g.\n',n/p, SparseLev(k)/num_groups);
            dopts.seednum = dopts.seednum + k*l;        % set seed number            
            % Generate data;                       
            dopts.matrixtype='logistic';          
            [A,At,b,be,xe,supp,suppg,gidx] = gendata(n,p,num_groups,num_nz_groups,dopts);
            % Create a group index cell array that accommodates both equal and unequal group sizes
            sgidx = arrayfun(@(kki) find(gidx == kki), (1:num_groups)', 'UniformOutput', false);
           
            
            %% HIWT-GSC
            ff = 'HIWT_GSC';           
            if ismember(ff,FIELDS)
                method_no = 1;
                x0=zeros(p,1);
                opts_HIWT.sgidx = sgidx;
                opts_HIWT.gidx = gidx;
                opts_HIWT.s = num_nz_groups; % the desired cardinality (i.e., the number of nonzero groups)
                opts_HIWT.num_groups = num_groups;
                opts_HIWT.app = 'LR';
                Flogistic = @(x)myLogistic(x,At,be);
                printf('\n-- %s, begin at %s --\n',ff,datestr(now));
                tic;
                result_HIWTGSC = HIWT_GSC(Flogistic,A,be,x0,opts_HIWT);
                Scputime(l,k,method_no) = toc;
                printf('-- %s is done, at %s --\n',ff,datestr(now));
                P = length(intersect(result_HIWTGSC.T,suppg))/length(result_HIWTGSC.T);
                R = length(intersect(result_HIWTGSC.T,suppg))/num_nz_groups;
                F1(l,k,method_no) = 2*P*R/(P+R);
                q = 1./(1+exp(-A*result_HIWTGSC.x));
                y=zeros(n,1);
                for i = 1:n    
                    if q(i)>0.5
                        y(i)=1;
                    else
                        y(i)=0;
                    end   
                end
                ER(l,k,method_no) = length(find((y-be)~=0))/n;  % error rate
            end
                   
        end       
end


    %% Compute average performance
    printf('\n=========================Average Results=========================\n');
    printf('n = %d, p = %d, Kg = %d, Sparse Level start at: %g, end at: %g\n',n, p, num_groups, SparseLev(1)/num_groups, SparseLev(len)/num_groups );

    average_time = mean(Scputime)
    average_F1 = mean(F1)
    average_ER = mean(ER)

 

