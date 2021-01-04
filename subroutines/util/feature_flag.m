function value = feature_flag(flag_name, default_value)
% function value = feature_flag(flag_name, default_value)
% 
% Return the value of the specified feature flag (environment variable), or the default value if
%   it is not set or empty
% %

value = getenv(flag_name);
if isempty(value)
    value = default_value;
end
