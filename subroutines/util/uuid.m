function x = uuid(varargin)
    % Return (leading characters of) a new UUID
    %
    % Arguments:
    %   n: Optional integer, number of (leading) hexadecimal characters of the UUID to return, if
    %       not specified you get the full 32, and if you ask for more than that you get an error
    %
    % see: https://www.briandalessandro.com/blog/how-to-create-a-random-uuid-in-matlab/
    % %
    
    assert(nargin <= 1, 'Expect at most one input argument, check function help for details');
   
    x = char(java.util.UUID.randomUUID);
    x = strrep(x, '-', '');  % Java library inserts hyphens for readability, we drop them
    
    
    if nargin == 1
        assert(varargin{1} <= 32, 'UUID has at most 32 hexadecimal charactors');
        x = x(1:varargin{1});
    end
        
end