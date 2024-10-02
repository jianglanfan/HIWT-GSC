% =========================================================================
% Project: HIWT-GSC
% Repository: https://github.com/<your-username>/HIWT-GSC
%
% File Name: HIWT_GSC.m
% Author: Lanfan Jiang
% Date Created: September 2024
% Last Modified: September 29, 2024
% Version: N/A
%
% Description:
%   This function implements the HIWT-GSC algorithm to solve 
%   group-sparsity-constrained optimization problems. It leverages a 
%   homotopy approach to gradually adjust the regularization parameter, 
%   ensuring better convergence properties when solving 
%   group-sparse-constrained optimization problems. 
%   During each iteration of the homotopy loop, the HIWT-GSC algorithm 
%   calls an inner optimization algorithm IWT-GSC to refine the solution 
%   at the current regularization level. It can be applied 
%   to various applications, such as least squares regression, logistic 
%   regression, and image classification.
%
% Usage:
%   homo_logger = HIWT_GSC(fun_obj, A, b, x0, options)
%
% Inputs:
%   fun_obj    - Objective function handle.
%   A          - Sample matrix.
%   b          - Sampled data.
%   x0         - Initial solution (defaults to 0).
%   options    - A structure containing:
%       .num_groups - Number of groups in the signal.
%       .s          - Estimated upper bound on the number of non-zero groups.
%       .num_stages - [default = 10] Number of homotopy iterations in HIWT-GSC.
%       .tol_x      - [default = 1e-3] Initial accuracy tolerance for IWT-GSC solutions.
%       .min_tolx   - [default = 1e-6] Minimal value for tol_x.
%       .de_tolx    - [default = 0.2] Decrease factor for tol_x.
%       .eta        - [default = 2] Increase factor for lambda.
%       .equalsize  - [default = true] Whether group sizes are equal.
%       .strategy   - [default = 'B'] Support set identification strategy:
%                     'B': Best-s Groups, 'T': Top-s Groups, 'H': Hybrid.
%       .app        - [default = 'LS'] Application type:
%                     'LS': Least Squares, 'LR': Logistic Regression,
%                     'IC': Image Classification.
%       .debias     - [default = true] Whether to perform a debias operation at the end.
%       .denoise    - [default = false] Whether to use 'residue < del' as a 
%                     halting condition for noisy data.
%       .gidx       - Group index formatted as:
%                     [ones(d_1,1); 2*ones(d_2,1); ... 
%                        num_groups*ones(d_{num_groups},1)] where d_i is 
%                        the length of group_i for i = 1,...,num_groups.
%
% Outputs:
%   x - Estimated solution.
%   T - Support set of x at the group level.
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

function homo_logger = HIWT_GSC(fun_obj,A,b,x0,options)

    % Validate mandatory fields
    required_fields = {'num_groups', 's', 'gidx'};
    validate_options(options, required_fields);  % Function for required field checks
    
    % Default values of the fields in 'options'
    flds = {'num_stages','de_tolx','eta','app','debias','denoise','min_tolx', 'num_groups', 's'};
    vals = { 10, 0.2, 2, 'LS', true, false, 1e-6, 0, 0}; 
    if exist('options','var')  
        for k = 1:numel(flds)
            if isfield(options,flds{k})
                vals{k} = options.(flds{k});
            end
        end
    end
    [ num_stages, de_tolx, eta, app, debias, denoise, min_tolx, num_groups, s] = deal(vals{:});
    
    
    if ~isfield(options,'lambda')
        options.lambda = 0.1; 
    end
    if ~isfield(options,'tol_x')
        options.tol_x = 1e-3; 
    end
    if ~isfield(options,'strategy')
        options.strategy = 'B'; 
    end
    if ~isfield(options,'tau')
        options.tau = 1; 
    end
    if ~isfield(options,'sgidx')
        sgidx = arrayfun(@(kki) find(options.gidx == kki), (1:num_groups)', 'UniformOutput', false);
        options.sgidx = sgidx;
    end
        
    switch app
        case 'IC'                   % Image Classification
            de_tolx = 1;
            denoise = true;
            options.strategy = 'H'; % Hybrid strategy
        case 'LR'
            debias = false;         % We do not provide debias code for logistic regression in this version
    end
    delta_s = ceil(s/4);
    options.s = delta_s;
    p = size(A,2);  
    options.p = p;
    x = x0;
    w = ones(p,1);
    
    for k = 1:num_stages                 
        logger = IWT_GSC(fun_obj,x,w,options);
        x = logger.x;
        w = logger.w;
        g = logger.g;
        options.tau = logger.tau;       

        % halting condition
        halt_cond = (options.s == s);
        halt_cond1 = (sum(abs(w.*x)) == 0);        
        residue = norm(g,inf)/options.lambda;
        halt_cond2 = (residue < 1e-3);  
        if (halt_cond && halt_cond1 && halt_cond2)
            break;
        end
        if denoise
            if num_groups ~= p
                ind = vertcat(options.sgidx{logger.T});              
            else
                ind = logger.T;
            end
            x_tmp = zeros(p,1);
            B = A(:,ind);
            x_tmp(ind) = (inv(B'*B))*B'*b;  
            residue = norm(B*x_tmp(ind)-b);              
            halt_cond3 = (residue<options.del);
            if halt_cond3
                break;
            end
        end
        
        options.lambda = eta*options.lambda; %increase lambda in each iteration
        % Pass different tolerance as a parameter of IWT_GSC in each iteration
        options.tol_x = max(min_tolx, options.tol_x*de_tolx);
        options.s = min(s,ceil(options.s*2));
    end
  
    %% debias
    if debias == true  
        if num_groups ~= p
            ind = vertcat(options.sgidx{logger.T});              
        else
            ind = logger.T;
        end
        switch app
            case 'LS'
                x = zeros(p,1);
                B = A(:,ind);
                x(ind) = (inv(B'*B))*B'*b; 
            case 'IC'
                x = zeros(p,1);
                B = A(:,ind);
                x(ind) = (inv(B'*B))*B'*b; 
            case 'LR'  %logistic
                disp('The authors do not provide debias code for logistic regression in this version.')
            otherwise
                disp('Undefined type')
        end
    end
    
    %% update solution
    homo_logger.x = x; 
    homo_logger.T = sort(logger.T);    
end



