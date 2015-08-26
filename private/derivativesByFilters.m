function [varargout] = derivativesByFilters(Fxy, varargin)
%derivativesByFilters - 1st and 2nd discrete derivatives from 5-tap and 7-tap filters
%
% This function computes 1st and 2nd derivatives of a scalar field F(x,y) using
% a differential filter with 5 or 7 tap coefficients. Farid and Simoncelli (2004)
% introduced an optimal filter design that gives results that are significantly 
% more accurate than MATLAB's GRADIENT function for gradients that are at angles
% other than vertical or horizontal. This in turn greatly improves estimates of 
% gradient orientation.
%
% [gx, gy, gxx, gyy, gxy] = derivativesByFilters(Fxy, derivative specifiers)
% The input arguments are:
%    Fxy: F(x,y) is the scalar field to be differentiated.
%    derivative specifiers: A comma-separated list of character strings
%       from the set {'x', 'y', 'xx', 'yy', 'xy'}, representing the
%       derivatives gx, gy, gxx, gyy, and gxy (as described below).
%       The specifiers can be in any order.
%    gx, gy, gxx, gyy, gxy: 1st derivatives in x and y, 2nd derivatives
%       in x and y, and the cross derivative relative to x and y, respectively.
%       The order of these output arguments will match the order of the
%       derivative specifiers in the input arguments.
%
% [..] = derivativesByFilters(..., h), where h specifies the grid size for
%    a uniform grid. The default is h = 1;
%
% [..] = derivativesByFilters(..., hx hy), where hx and hy specify grid size 
%    in the x and y directions, respectively. The default is hx = hy = 1.
%
% [..] = derivativesByFilters(..., filterSize), where filterSize is a string
%    that specifies the size of the filter. The default is 'five', and
%    alternative is 'seven', refering to 5-tap and 7-tap filters, respectively.
%    The 7-tap filter is more precise, but otherwise requires more cpu time
%    and produces estimates with lower spatial resolution.
%
% Examples:
%   Compute 1st derivatives in x and y:
%   [gx, gy] = derivativesByFilters(Fxy, 'x', 'y');
%
%   Compute 2nd derivative in x, and 1st derivative in y:
%    [gxx, gy] = derivative5(Fxy, 'xx', 'y')
%
% Note that derivative grids are returned with the same size as
% the input grid Fxy. The function uses Matlab's conv2 function, and,
% as a result, the convolution of the filters with Fxy are done using
% zero-padding. Other kinds of padding can be used by applying the
% padding to Fxy before the call, and removing the effects after
% the call. Padding for the 5-tap and 7-tap filter options
% requires the addition of 2 and 3 elements, respectively, to
% the perimeter of Fxy.
%
% See also: GRADIENT
%
% Reference: Hany Farid and Eero Simoncelli "Differentiation of Discrete
% Multi-Dimensional Signals" IEEE Trans. Image Processing. 13(4): 496-508 (2004)

% Copyright (c) 2010 Peter Kovesi
% Centre for Exploration Targeting
% The University of Western Australia
% http://www.csse.uwa.edu.au/~pk/research/matlabfns/
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.
%
% April 2010, Peter Kovesi (original code is available at web link above)
% August 2012, Mark Brandon: minor modifications, including arguments for
% hx and hy, and consolidation of the original functions derivative5.m
% and derivative7.m into a single function, derivativesByFilter.m.

%% Start of function
% Parse input and output variables
% Check input arguments for grid-size variables hx, or hy, and
% filterSize.
N = length(varargin);
hx = 1;
hy = 1;
filterSize = 'five';
if N>0 && iscellstr(varargin(N))
    %... Set to selected filterSize
    filterSize = char(varargin(N));
    N = N-1;
    if ~strcmp(filterSize,'five') && ~strcmp(filterSize,'seven')
        error('derivativesByFilter: incorrect argument, ''%s'', for filterSize',filterSize);
    end
end
if isnumeric(varargin{N})
    hy = varargin{N};
    N = N - 1;
    if isnumeric(varargin{N})
        hx = varargin{N};
        N = N - 1;
    else
        hx = hy;
    end
end
%... Trim optional arguments, but leave derivative specifiers
varargout = cell(size(1,N));

% Check if we are just computing 1st derivatives.  If so use the
% interpolant and derivative filters optimized for 1st derivatives, else
% use 2nd derivative filters and interpolant coefficients.
% Detection is done by seeing if any of the derivative specifier
% arguments is longer than 1 char, this implies 2nd derivative needed.
calc2ndDerivative = false;
for i = 1:N
    if length(varargin{i}) > 1
        calc2ndDerivative = true;
        break
    end
end

%... Construct tap filters
switch filterSize
    case 'five'
        if ~calc2ndDerivative
            % 5-tap 1st-derivative cofficients. These are optimal if you are just
            % seeking the 1st deriavtives
            p = [0.037659  0.249153  0.426375  0.249153  0.037659];
            d1 =[0.109604  0.276691  0.000000 -0.276691 -0.109604];
        else
            % 5-tap 2nd-derivative coefficients. The associated 1st-derivative
            % coefficients are not quite as optimal as the ones above but are
            % consistent with the 2nd-derivative interpolator p and thus are
            % appropriate to use if you are after both 1st and 2nd derivatives.
            p  = [0.030320  0.249724  0.439911  0.249724  0.030320];
            d1 = [0.104550  0.292315  0.000000 -0.292315 -0.104550];
            d2 = [0.232905  0.002668 -0.471147  0.002668  0.232905];
        end       
    case 'seven'
        % 7-tap interpolant and 1st- and 2ndderivative coefficients.
        p  = [ 0.004711  0.069321  0.245410  0.361117  0.245410  0.069321  0.004711];
        d1 = [ 0.018708  0.125376  0.193091  0.000000 -0.193091 -0.125376 -0.018708];
        d2 = [ 0.055336  0.137778 -0.056554 -0.273118 -0.056554  0.137778  0.055336];
    otherwise
        error('derivativeByFilters: incorrect argument, ''%s'', for filterSize.',filterSize)
end

% Compute derivatives. Note that the convolution is calculated in conv2 
% in two step, with a 1D convolution down the columns, followed by a 1D
% convolution along the rows.
dFdxExisits = false;
for i = 1:N
    switch logical(true)
        case strcmpi('x', varargin{i})
            %... Calculate dF/dx
            varargout{i} = conv2(p, d1, Fxy, 'same')/hx;
            % Set flag dfdxExisit to true for calculation of dF/dxdy
            dFdxExisits = true;
            idFdx = i;
        case strcmpi('y', varargin{i})
            %... Calculate dF/dy
            varargout{i} = conv2(d1, p, Fxy, 'same')/hy;
        case strcmpi('xx', varargin{i})
            %... Calculate d2F/dx2
            varargout{i} = conv2(p, d2, Fxy, 'same')/hx^2;
        case strcmpi('yy', varargin{i})
            %... Calculate d2F/dy2
            varargout{i} = conv2(d2, p, Fxy, 'same')/hy^2;
        case strcmpi('xy', varargin{i}) || strcmpi('yx', varargin{i})
            %... Calculate d2F/dxdy
            if dFdxExisits
                varargout{i} = conv2(d1, p, varargout{idFdx}, 'same')/hy;
            else
                dFdx = conv2(p, d1, Fxy, 'same')/hx;
                varargout{i} = conv2(d1, p, dFdx, 'same')/hy;
            end
        otherwise
            error('''%s'' is an unrecognized derivative option',varargin{i});
    end
end
end
