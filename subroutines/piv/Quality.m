classdef Quality < uint8
    % Enumeration specifying PIV quality flags
    
    enumeration
        Valid (0)
        EmptySampleWindow (1)
        BelowMinOverlap (2)
        ValidationFailed (3)
        PeakFindingFailed (4)
        Skipped (5)
    end
    
    methods(Static)
        
        function [uint8_array] = to_uint8(obj_array)
            % cast array of Quality objects to array of uint8 values
            uint8_array = uint8(obj_array);
        end
        
    end
    
end     
    
    
        