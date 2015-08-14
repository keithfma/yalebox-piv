function [gx gy] = nangrad2(f,h)
% function [gx gy] = nangrad2(f,h)
%
% Numerical gradient ignoring NaNs.  Computes central differences (2nd
% order acc) where possible and forward or backward differences (1st order
% acc) where not possible.  If input argument h is not provided, unit grid
% spacing is assumed.. Useful to compute gradients in irregular regions
% (bounded by NaNs).

if nargin<2
    h = 1;
end

[nr nc] = size(f);
gx = nan(nr,nc);
gy = nan(nr,nc);

% compute x-direction gradient (along rows)
for r = 1:nr
    
    % left endpoint, forward difference
    if ~isnan(f(r,2)) && ~isnan(f(r,1)) % data exists, compute gradient
        gx(r,1) = f(r,2)-f(r,1);
    else
        gx(r,1) = NaN;
    end
    
    % interior points, depending on data availability
    for c = 2:nc-1
       
        if ~isnan(f(r,c)) && (~isnan(f(r,c-1)) || ~isnan(f(r,c+1))) % data exists compute gradient 
            
            if ~isnan(f(r,c-1)) && ~isnan(f(r,c+1)) % central difference
               gx(r,c) = ( f(r,c+1)-f(r,c-1) )/2;
               
            elseif ~isnan(f(r,c-1)) && isnan(f(r,c+1)) % backward difference
                gx(r,c) = f(r,c)-f(r,c-1);
               
            elseif isnan(f(r,c-1)) && ~isnan(f(r,c+1)) % forward difference
                gx(r,c) = f(r,c+1)-f(r,c);
            end
        
            
        else
            gx(r,c) = NaN;
        end
        
        
    end
    
    % right endpoint, backward difference
    if ~isnan(f(r,nc)) && ~isnan(f(r,nc-1)) % data exists, compute gradient
        gx(r,nc) = f(r,nc)-f(r,nc-1);
    else
        gx(r,nc) = NaN;
    end
end


% compute y-direction gradient (along columns)
for c = 1:nc
    
    % lower endpoint, forward difference
    if ~isnan(f(2,c)) && ~isnan(f(1,c)) % data exists, compute gradient
        gy(1,c) = f(2,c)-f(1,c);
    else
        gy(1,c) = NaN;
    end
    
    % interior points, depending on data availability
    for r = 2:nr-1
       
        if ~isnan(f(r,c)) && (~isnan(f(r-1,c)) || ~isnan(f(r+1,c))) % data exists compute gradient 
            
            if ~isnan(f(r-1,c)) && ~isnan(f(r+1,c)) % central difference
               gy(r,c) = ( f(r+1,c)-f(r-1,c) )/2;
               
            elseif ~isnan(f(r-1,c)) && isnan(f(r+1,c)) % backward difference
                gy(r,c) = f(r,c)-f(r-1,c);
               
            elseif isnan(f(r-1,c)) && ~isnan(f(r+1,c)) % forward difference
                gy(r,c) = f(r+1,c)-f(r,c);
            end
        
            
        else
            gy(r,c) = NaN;
        end
        
        
    end
    
    % right endpoint, backward difference
    if ~isnan(f(nr,c)) && ~isnan(f(nr-1,c)) % data exists, compute gradient
        gy(nr,c) = f(nr,c)-f(nr-1,c);
    else
        gy(nr,c) = NaN;
    end
end

% adjust spacing
gx = gx/h;
gy = gy/h;

