% =========================================================================
% Project: HIWT-GSC
% Repository: https://github.com/jianglanfan/HIWT-GSC
%
% File Name: gendata.m
% Author: Bangti Jin (bangti.jin@gmail.com), 
%         Yuling Jiao (yulingjiaomath@whu.edu.cn)
% Date Created: September 4, 2015
% Last Modified: September 2024 by Lanfan Jiang
% Version: N/A
%
% Description:
%   This function generates one-dimensional simulated data based on the 
%   specified parameters, including the number of samples, signal length, 
%   and group characteristics.
%
% Usage:
%   [A, At, b, be, x, supp, suppg, gidx] = gendata(n, p, num_groups, 
%   num_nz_groups, opts)
%
% Inputs:
% - n             - number of samples
% - p             - signal length
% - num_groups    - number of groups in the signal
% - num_nz_groups - number of nonzero groups in the signal
% - opts          - structure containing the following fields:
%     .ratio      - range of value in the signal (= 10^ratio) 
%                    (default: 1)
%     .sigma      - noise variance (default: 0.1)
%     .seednum    - seed number (default: 0)
%     .gsize      - integer vector of length num_groups containing sizes 
%                   of g_i (default: [], for equal group size)
%     .isnorm     - normalization after adding correlation 
%                   (default: 1)
%     .matrixtype - type of sample matrix ('gauss' (default) or 
%                   'logistic')
%
% Outputs:
% - A            - normalized sample matrix
% - At           - transpose of A
% - b            - sampled data with noise
% - be           - sampled data without noise
% - x            - the generated signal
% - supp         - support of the generated signal x
% - suppg        - support of x in group level
% - gidx         - group index formed as
%                  [ones(d_1,1); 2*ones(d_2,1); ...; 
%                  num_groups*ones(d_{num_groups},1)]
%                  with d_i being the length of group_i, i=1,...,num_groups
%
% =========================================================================


function [A,At,b,be,x,supp,suppg,gidx] = gendata(n,p,num_groups,num_nz_groups,opts)

disp(' Generating data ...')

% Initialize parameters
if isfield(opts,'sigma')
    sigma = opts.sigma;
else
    sigma = .1; 
end
if isfield(opts,'ratio')
    ratio = opts.ratio;
else
    ratio = 1;
end
if isfield(opts,'seednum')
    seednum = opts.seednum;
else
    seednum = 0; 
end
if isfield(opts,'matrixtype')
    matrixtype = opts.matrixtype;
else
    matrixtype = 'gauss';
end
if isfield(opts,'isnorm')
    isnorm = opts.isnorm;
else
    isnorm = 1;
end
if isfield(opts,'gsize')
    gsize = opts.gsize;
else
    gsize = [];
end
rng(seednum);

% generate signal  
x = zeros(p,1);                     % true signal
if isempty(gsize) 
    avgsize = round(p/num_groups);	% average group size make
    gidx = [];                      % index for groups
    for k = 1:num_groups  
        gidx = [gidx; k*ones(avgsize,1)];   
    end
    gidx = [gidx;num_groups*ones(p - avgsize*num_groups,1)];  % in case p/num_groups is not an integer
else
    if sum(gsize) ~= p
        disp('error, the sum of goup size must equal to p')
    else
        gidx = [];              
        for k = 1:num_groups
            gidx = [gidx;k*ones(gsize(k),1)];   
        end
    end
end
q  = randperm(num_groups);              % random permutation of num_groups
suppg = sort((q(1:num_nz_groups)))';    % group support 
idx = ismember(gidx,suppg);
supp = find(idx);           % find the support of xe
K = length(supp);           % return the number of nonzero elements of xe 

if ratio ~= 0
    vx = ratio*rand(K,1);
    vx = vx-min(vx);
    vx = vx/max(vx)*ratio;
    x(supp) = 10.^vx.*sign(randn(K,1));
else
    for k = 1:length(suppg)
        id = find(gidx == suppg(k));
        x(id) = k*sign(randn(length(id),1));
    end
end
    

% generate matrix A
switch matrixtype
    case 'gauss'
        A = randn(n,p);
        A = normalize(A,isnorm);
        At = A';
        be  = A*x;
        b = be+sigma*randn(n,1);
    case 'logistic'
        A = randn(n,p);
        A = normalize(A,isnorm);
        At = A';    
        q = 1./(1+exp(-A*x));   % intercept = 0 
        y = zeros(n,1);
        for i = 1:n    
            if q(i)>0.5
                y(i)=1;
            else
                y(i)=0;
            end
        end
        be = y;
        b = be;
    otherwise
        disp('Undefined matrix type')
end

disp(' Data generation is done!')
end


function [sX] = normalize(X,isnorm)
%------------------------------------------------------------------------%
% purpose: normalize the sensing matrix X  %(columnwise)                 %                                                  %
% Inputs:                                                                %
%    X     --- un-normalized matrix                                      %                                    %
%  isnorm  --- 1, if normalized                                          %
% Outputs:                                                               %
%     sX     --- normalized matrix                                       %      
%------------------------------------------------------------------------%
p  = size(X,2);
sX = X;

% normalization,sX(:,k)'*sX(:,k)=1
if isnorm == 1  
    for k =1:p
        sX(:,k) = sX(:,k)/norm(sX(:,k));
    end
end
end


