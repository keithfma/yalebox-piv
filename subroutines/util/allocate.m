function allocate(file_obj, var_name, data_type, dimensions)
% function allocate(file_obj, var_name, data_type, dimensions)
%
% Allocate space in matfile output variable
% 
% Arguments:
%   file_obj: matfile object, open for writing
%   data_type: function handle for the variable data type, e.g., @double
%   var_name: string, name of variable to allocate in matfile
%   dimensions: 1D array, size of variable to allocate
% %

switch data_type
    
    case 'double'
        empty = @double.empty;
        last = NaN;

    case 'single'
        empty = @single.empty;
        last = single(NaN);
        
    case 'uint8'
        empty = @uint8.empty;
        last = uint8(0);
    
    case 'logical'
        empty = @logical.empty;
        last = false;
        
    otherwise
        error('Bad value for data_type: %s', data_type);
end

file_obj.(var_name) = empty(zeros(size(dimensions)));
dimensions = num2cell(dimensions);
file_obj.(var_name)(dimensions{:}) = last;
