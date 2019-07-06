classdef Quality < uint8
    % Enumeration specifying PIV quality flags
    
    enumeration
        Valid (0)
        EmptySampleWindow (1)
        BelowMinOverlap (2)
        ValidationFailed (3)
    end
    
    methods(Static)
        
        function [valid_array] = is_valid(obj_array)
            % return true where input array is Quality.Valid, false elsewhere
            valid_array = obj_array == Quality.Valid;
        end
        
        function [uint8_array] = to_uint8(obj_array)
            % cast array of Quality objects to array of uint8 values
            uint8_array = uint8(obj_array);
        end
        
    end
    
end     
    
    
        