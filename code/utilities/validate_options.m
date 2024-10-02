% =========================================================================
% Project: HIWT-GSC
% Repository: https://github.com/jianglanfan/HIWT-GSC
%
% File Name: validate_options.m
% Author: Lanfan Jiang
% Date Created: September 2024
% Last Modified: September 30, 2024
% Version: N/A
%
% Description:
%   This function checks whether all the required fields are present in the 
%   user-provided options structure. If any required field is missing, it 
%   throws an error indicating the missing field.

% Usage:
%   validate_options(options, required_fields)
%   - options: Structure containing the parameters passed to the solver.
%   - required_fields: Cell array of strings, listing the required fields 
%                      that must be present in the options structure.
%
% Inputs:
% - options: Structure containing user-specified parameters.
% - required_fields: Cell array listing required fields as strings.
%
% Outputs:
% - None. Throws an error if any required field is missing.
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

function validate_options(options, required_fields)
    for i = 1:numel(required_fields)
        if ~isfield(options, required_fields{i})
            error('Missing required field: %s', required_fields{i});
        end
    end
end