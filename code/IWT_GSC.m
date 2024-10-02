% =========================================================================
% Project: HIWT-GSC
% Repository: https://github.com/jianglanfan/HIWT-GSC
%
% File Name: IWT_GSC.m
% Author: Lanfan Jiang
% Date Created: September 2024
% Last Modified: September 29, 2024
% Version: N/A
%
% Description:
%   This function implements the Iterative-Weighted Thresholding (IWT) 
%   algorithm for solving group-sparsity-constrained optimization problems.
%
% Usage:
%   [x, T, w, g, tau] = IWT_GSC(fun_obj, x, w, options)
%
% Inputs:
% - fun_obj    - Objective function handle
% - x          - Initial solution
% - w          - Initial weight
% - options    - A structure variable determining various options of
%                the algorithm through the following fields:
%       .strategy       - [default = 'B']: Support set identification 
%                         strategy. Takes values 'B', 'T', and 'H':
%                         'B': Best-s groups
%                         'T': Top-s groups
%                         'H': Hybrid strategy
%                         You can customize different support set identification 
%                         strategies for different applications.
%       .s              - An integer determining the desired sparsity level
%       .x_norm         - [default = 2]: 2 for L_{2,1}-norm, 1 for 
%                         L_{1,1}-norm
%       .tau            - [default = 1]: Initial stepsize
%       .stepsizeShrink - [default = 0.5]: Factor by which stepsize tau 
%                         is reduced: tau = stepsizeShrink * tau
%       .maxIter        - [default = 200]: The maximum number of the 
%                         IWT_GSC iterations to be performed
%       .tol_x          - [default = 1e-3]: IWT_GSC halts if the 
%                         decrease in the value of the relative error 
%                         is less than tol_x
%       .gamma          - [default = 0.7]: Controls the degree of 
%                         non-monotonicity
%
% Outputs:
% - x          - Estimated solution
% - T          - Support of x in group level
% - w          - Weight; the selected groups in T are weighted with 0
% - g          - Gradient
% - tau        - Stepsize
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

function logger = IWT_GSC(fun_obj,x,w,options)

    p = options.p;  % p = size(x,1);
    % Default values of the fields in 'options'
    flds = { 's','sgidx', 'gidx', 'num_groups', 'app','lambda','tau','taumax','taumin','stepsizeShrink','tol_x','maxIter','equalsize','x_norm','strategy','verbose'}; 
    vals = { 0, [], [], p, 'LS', 1e-1, 1, 1e+8, 1e-8, 0.5, 1e-3, 200, true, 2, 'B', false};   
    if exist('options','var')  
        for k = 1:numel(flds)
            if isfield(options,flds{k})
                vals{k} = options.(flds{k});
            end
        end
    end
    [s, sgidx, gidx, num_groups, app, lambda, tau, taumax, taumin, stepsizeShrink, tol_x, maxIter, equalsize, x_norm, strategy, verbose] = deal(vals{:}); % Set the options of IWT_GSC 
     
    % gamma is used to tune non-monotonicity during line searches
    if isfield(options,'gamma')
        gamma = options.gamma;
    else
        gamma=0.7;  
    end
    
      
    %% Initialization
    if equalsize
        gs = p/num_groups;	% group size: number of elements in each group
    end
    T = [];                 % support set 
	[fx, g] = fun_obj(x);
    Rx = get_Rx(x, w, num_groups, p, x_norm, gidx, sgidx, equalsize);     
    Fx = fx+lambda*Rx;
    ck = Fx; thetak = 1; sigma = 5e-10;

        
    %% Main Iterations of the IWT_GSC Algorithm
    for j = 1:maxIter      % inner iteration number(exclude line search times)
        % Update the Record
        Fx_old = Fx;
        x_old = x; 
        g_old = g;
        tau_old = tau;
        T_old = T;
        
        % line search  
        backtrackCount = 0; % maximum number of backtracking loops for line searches
        ls_pass = 0;        % sign of a successful line search
        grouped_z = zeros(num_groups,1);
        while (tau>taumin) && (backtrackCount<20)  % backtracking loop
            backtrackCount = backtrackCount+1;
            z = x_old - tau*g_old; 
            z_abs = abs(z);
            nu = lambda*tau;    % threshold value
             % compute soft(z)
            if x_norm == 2
                z_square = z.^2;
                grouped_z = grouped_value(z_square, gidx, sgidx, equalsize);
                zg_norm = sqrt(grouped_z);      %||z_G||_{2,1}
                gc = max(0,zg_norm-nu);
                gc = gc./max(realmin,zg_norm);  % group coefficient=(zg_norm-nu)/zg_norm
                if equalsize
                    gca = repelem(gc,gs);       % group coefficient append
                else
                    len = cellfun(@length, sgidx);  % get the length of each sgidx
                    gca = zeros(p, 1);  
                    start_idx = 1;
                    for kki = 1:num_groups
                        end_idx = start_idx + len(kki) - 1;
                        gca(start_idx:end_idx) = gc(kki);  
                        start_idx = end_idx + 1;
                    end
                end
                soft_z = gca.*z;               
            else
                if x_norm == 1                   
                    soft_z = sign(z).*max(0,z_abs-nu);	
                end
            end
    
            % Identify support set, you can add your own strategy here
            switch strategy                   
                case 'B' % Best-s groups
                    if num_groups ~= p
                    	if x_norm == 2 
                           lx = 0.5*((soft_z-z).^2)/tau;
                           soft_z_square = soft_z.^2;
                           grouped_lx = grouped_value(lx, gidx, sgidx, equalsize);
                           grouped_rx = grouped_value(soft_z_square, gidx, sgidx, equalsize);
                           grouped_rx = sqrt(grouped_rx);
                           grouped_loss = grouped_lx+lambda*grouped_rx;
                        else 
                        	if x_norm == 1
                                loss_function_value = 0.5*((soft_z-z).^2)/tau+lambda*abs(soft_z);
                                grouped_loss = grouped_value(loss_function_value, gidx, sgidx, equalsize);  
                        	end
                    	end
                    else
                        grouped_loss = 0.5*((soft_z-z).^2)/tau+lambda*abs(soft_z);
                    end  
                    [~,sorted_loss_idx] = sort(grouped_loss,'descend');
                    T = sorted_loss_idx(1:min(numel(find(grouped_loss>0)),s));

                case 'T'  % Top-s groups
                    if num_groups ~= p
                    	if x_norm == 2
                            grouped_z = zg_norm;
                        else
                            if x_norm == 1
                                grouped_z = grouped_value(z_abs, gidx, sgidx, equalsize);
                            end
                    	end                       
                    else
                        grouped_z = z_abs;
                    end 
                    [~,sorted_groupedz_idx] = sort(grouped_z,'descend');
                    T = sorted_groupedz_idx(1:min(numel(find(grouped_z>0)),s));

                case 'H'  % Hybrid strategy
                    if num_groups ~= p
                        if x_norm == 2                         
                           lx = 0.5*((soft_z-z).^2)/tau;
                           soft_z_square = soft_z.^2;
                           grouped_lx = grouped_value(lx, gidx, sgidx, equalsize);
                           grouped_rx = grouped_value(soft_z_square, gidx, sgidx, equalsize);
                           grouped_rx = sqrt(grouped_rx);
                           grouped_loss = grouped_lx+lambda*grouped_rx;
                        else 
                            if x_norm == 1
                                loss_function_value = 0.5*((soft_z-z).^2)/tau+lambda*abs(soft_z);
                                grouped_loss = grouped_value(loss_function_value, gidx, sgidx, equalsize);  
                            end
                        end
                    else
                        grouped_loss = 0.5*((soft_z-z).^2)/tau+lambda*abs(soft_z);
                    end                     
                    [~,sorted_loss_idx] = sort(grouped_loss,'descend'); 
                    [~,sorted_z_idx] = sort(z_abs,'descend');
                    T1 = sorted_loss_idx(1:min(numel(find(grouped_loss>0)),s));
                    T2_tmp = sorted_z_idx(1:min(numel(find(z_abs>0)),s));
                    T2 = gidx(T2_tmp);
                    T_tmp = [T1' T2'];
                    T = unique(T_tmp);   
                otherwise
                    disp('Undefined strategy')
                    return;
            end

            % determin x_Gi and w_i according to group support set T
            w = ones(p,1);
            if s ~= 0
                if num_groups ~= p
                    ind = vertcat(sgidx{T});              
                else
                    ind = T;
                end
            	x = soft_z;
            	w(ind) = 0; 
            	x(ind) = z(ind);
            else
                x = zeros(p,1);
            end
           
          %% updata fx,g,Fx,ck
            [fx, g] = fun_obj(x);
            Rx = get_Rx(x, w, num_groups, p, x_norm, gidx, sgidx, equalsize);
            Fx = fx+lambda*Rx;
            Fx_threshold = ck-(sigma*(norm(x-x_old)^2)/tau);
            if Fx <= Fx_threshold
               ls_pass = 1;         % find an aprropriate stepsize
               break;
            else
               tau = max(taumin,tau*stepsizeShrink);
               if verbose
                    fprintf(1,'inner interation = %4d, backtrackCount = %2d, fx = %10.6e, decreasing stepsize to %6.2e\n', j, backtrackCount, fx, tau);
               end
            end          
        end
        
        
        if ~ls_pass
            x = x_old;
            g = g_old;
            Fx = Fx_old;
            T = T_old;
            tau = tau_old;
        end
        
        %% Check if Halting Condition Holds   
        switch app
            case 'LR'
            	residue =  abs(Fx_old - Fx); 
            	HaltCond = (residue<tol_x);
            otherwise
             	relerr = norm(x-x_old)/norm(x_old+realmin);
            	HaltCond = (relerr<tol_x);
        end        
        if HaltCond
            break;
        end
                
        % update ck
        thetak = 1+gamma*thetak;
        ck = ((thetak-1)*ck+Fx)/thetak;
        if verbose 
            fprintf(1,'find an aprropriate stepsize, inner interation=%4d, obj=%10.6e, stepsize=%e\n', j, fx, tau);  
        end
        
        %% BB
        dx = x-x_old;
        dg = g-g_old;
        dotprod = real(dot(dx,dg));            
        tau = (norm(dx)^2)/ (dotprod+realmin);  
        if tau <=0 || isinf(tau) || isnan(tau)     %  Make sure step is non-negative
            	tau = tau_old*1.5;  % let tau grow, backtracking will kick in if stepsize is too big
        end
        
        tau = max(taumin,min(taumax,tau));  
        
    end

    logger.x = x;
    logger.w = w;
    logger.T = T;
    logger.g = g;
    logger.tau = tau;     
end

function Rx_value = get_Rx(x, w, num_groups, p, x_norm, gidx, sgidx, equalsize) 
    if p ~= num_groups 
        switch x_norm
            case 1
                Rx_value = sum(abs(w.*x));
            case 2
                tmp = w.*x;
                x_square = tmp.^2;   
                if equalsize
                    grouped_x_square = accumarray(gidx,x_square);
                else
                    grouped_x_square = cellfun(@(ind) sum(x_square(ind)), sgidx, 'UniformOutput', true);  
                end
                Rx_value = sum(sqrt(grouped_x_square));
            otherwise
                disp('Illegal forms of norm')
                return;
        end
    else
        Rx_value = sum(abs(w.*x)); 
    end
end

% calculate the sum of the 'value' in group level
function g_value= grouped_value(value, gidx, sgidx, equalsize) 
    if equalsize
        g_value = accumarray(gidx, value);
    else
        g_value = cellfun(@(ind) sum(value(ind)), sgidx, 'UniformOutput', true); 
    end
end