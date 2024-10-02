% =========================================================================
% Project: HIWT-GSC
% Repository: https://github.com/jianglanfan/HIWT-GSC
%
% File Name: demo_LeastSquares.m
% Author: Lanfan Jiang
% Date Created: September 2024
% Last Modified: September 29, 2024
% Version: N/A
%
% Description:
%   This experiment demo shows how to use the HIWT-GSC for a 
%   group-sparsity-constrained least squares regression problem.
%   We evaluate the performance of HIWT-GSC in terms of probability 
%   of exact support recovery (PSR), relative error, and CPU time.
%
% Usage:
%   To run the demo, simply execute the corresponding script in MATLAB:
%   >> demo_LeastSquares
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
    p = 800;                % signal length  
    n = 200;                % number of samples, sampling rate = 25%
    num_groups = 200;       % number of groups in the signal 
    SparseLev = 20:2:36;    % the number of nonzero groups, sparse level: 10%-18%
else
    p = 10000;              % signal length
    n = 2000;               % number of samples, sampling rate = 20%
    num_groups = 2500;      % number of groups in the signal 
    SparseLev = 200:25:400; % the number of nonzero groups, sparse level: 8%-16%
end
gs = p/num_groups;          % number of elements in each group
n_method = numel(FIELDS);   % number of methods involved in testing
maxnum_test = 100;          % number of tests performed for each parameter setting, you can set maxnumtest = 10 for a quick test
len = length(SparseLev);

%% Define arrays to save the results of each test
Scputime = zeros(maxnum_test,len,n_method);
Srel2error = zeros(maxnum_test,len,n_method);   % relative error
PSR = zeros(maxnum_test,len,n_method);    % probability of exact support recovery, measured by S_G(x_true) = S_G(x_estimated).

dopts.sigma = 1e-1;         % noise variance (default: .1)
dopts.seednum  = 0;         % seed number   
rng('default'); 
noise = true;               % whether to use noisy data
fig = 1;
fid = 1;
printf = @(varargin) fprintf(fid,varargin{:});

for k = 1:len
    num_nz_groups = SparseLev(k);
    for l = 1:maxnum_test
        printf('\nSampling Rate = : %g, Sparse Level = %g.\n',n/p, SparseLev(k)/num_groups);
        dopts.seednum = dopts.seednum + k*l;               
        % Generate data;                       
        dopts.matrixtype='gauss';          
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
            opts_HIWT.s = num_nz_groups;    % the desired cardinality (i.e., the number of nonzero groups)
            opts_HIWT.num_groups = num_groups;
            opts_HIWT.gamma = 0.9;          % gamma is used to tune non-monotonicity during line searches
            printf('\n-- %s, begin at %s --\n',ff,datestr(now));
            tic;
            if noise == true
                Flinear = @(x)myLinear(x,A,b);
                opts_HIWT.denoise = true;
                opts_HIWT.del = norm(b-be);
                result_HIWTGSC = HIWT_GSC(Flinear,A,b,x0,opts_HIWT);
            else
                Flinear = @(x)myLinear(x,A,be);
                result_HIWTGSC = HIWT_GSC(Flinear,A,be,x0,opts_HIWT);
            end
            Scputime(l,k,method_no) = toc;
            printf('-- %s is done, at %s --\n',ff,datestr(now));
            Srel2error(l,k,method_no) = norm(result_HIWTGSC.x - xe)/norm(xe);
            if length(result_HIWTGSC.T) == length(suppg)
                if sort(result_HIWTGSC.T) == suppg
                    PSR(l,k,method_no) = 1;
                end
            end
        end

    end

    %% _______________________________  polt figs _______________________________
    % origianl signal
    figure(fig)
    fig = fig+1;
    scrsz = get(0,'ScreenSize');
    set(gcf,'Position',[10 10 0.9*scrsz(3) 0.9*scrsz(4)])
    subplot(2,1,1)
    plot(xe,'LineWidth',1.1)
    top = max(xe(:));
    bottom = min(xe(:));
    v = [0 p+1 bottom-0.05*(top-bottom)  top+0.05*((top-bottom))];
    set(gca,'FontName','Times')
    set(gca,'FontSize',14)
    if noise
        title(sprintf('Original (number of samples = %g, signal length = %g, number of groups = %g, active groups = %g, noise level = %g)',n,p,num_groups,num_nz_groups,dopts.sigma))
    else
        title(sprintf('Original (number of samples = %g, signal length = %g, number of groups = %g, active groups = %g)',n,p,num_groups,num_nz_groups))
    end
    axis(v)

    % reconstructed signal of HIWT-GSC
    subplot(2,1,2)
    plot(result_HIWTGSC.x,'LineWidth',1.1)
    top = max(result_HIWTGSC.x(:));
    bottom = min(result_HIWTGSC.x(:));
    v = [0 p+1 bottom-0.05*(top-bottom)  top+0.05*((top-bottom))];
    set(gca,'FontName','Times')
    set(gca,'FontSize',14)
    title(sprintf('result of HIWT-GSC: relative  error = %4e',Srel2error(l,k,1)))
    axis(v)
   %_____________________________ end of poltting figs _____________________________
        
end


 %% Compute average performance
printf('\n=========================Average Results=========================\n');
printf('n = %d, p = %d, Kg = %d, Sparse Level start at: %g, end at: %g\n',n, p, num_groups, SparseLev(1)/num_groups, SparseLev(len)/num_groups );
if noise
    printf('noise level = %g \n', dopts.sigma); 
end
average_time = mean(Scputime)
average_rel2error = mean(Srel2error)
average_PSR = mean(PSR)        


 

