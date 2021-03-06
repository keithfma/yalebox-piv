function [rimg, rx, ry] = movie_frame_resize(img, x, y, max_dim)
%
% Resample image such that it does not exceed the specified dimensions and the
% aspect ratio is preserved. Uses the MATLAB function imresize to interpolate
% and apply anti-aliasing filter. Requires that both output dimensions be even
% to be compatible with common video codecs (i.e. H.264)
%
% Arguments:
% 
%   rimg = 2D or 3D matrix, resized image 
%
%   rx, ry = Vectors, double, resized coordinate vectors corresponding to the
%       resized image.
%
%   img = 2D or 3D matrix, input image to be resized
%
%   x, y, = Vectors, double, coordinate vectors corresponding to the input
%       image
%
%   max_dim = 2-element vector, maximum dimensions for the resized output
%       image as [max_num_cols, max_num_rows]. 
%
% % Keith Ma, August 2015

nx = numel(x);
ny = numel(y);
nc = size(img,3); % colors

if nx > max_dim(1) || ny > max_dim(2) % resize
    
    [scale, ii] = min(max_dim(:)./[nx; ny]);
    if scale<1
        
        % resize to max width
        if ii==1 
            rimg = imresize(img, [NaN, max_dim(1)],...
                'Method', 'lanczos3',...            
                'Antialiasing', true);
        
        % resize to max height
        else 
            rimg = imresize(img, [max_dim(2), NaN],...
                'Method', 'lanczos3',...            
                'Antialiasing', true);
        end
        
        % resize coordinate vectors to match        
        rx = imresize(x(:), [size(rimg,2), 1],...
            'Method', 'lanczos3',...            
            'Antialiasing', true);            
        ry = imresize(y(:), [size(rimg,1), 1],...
            'Method', 'lanczos3',...            
            'Antialiasing', true);
        
    end
    
    % delete one row/col if dimension is not even
    if mod(numel(rx), 2) == 1
        rx = rx(1:end-1);
        rimg = rimg(:, 1:end-1, :);
    end
     if mod(numel(ry), 2) == 1
        ry = ry(1:end-1);
        rimg = rimg(1:end-1, :, :);
    end
    
    % trim values interpolated outside of valid [0, 1] range
    if isa(rimg, 'double')
        rimg(rimg>1) = 1;
        rimg(rimg<0) = 0;
    end
    
else % no resize needed
    rimg = img;
    rx = x;
    ry = y;
 
end
