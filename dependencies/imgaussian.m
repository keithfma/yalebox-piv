function I=imgaussian(I,sigma,filterSize)
% IMGAUSSIAN filters an 1D, 2D color/greyscale or 3D image with an
% Gaussian filter. This function is based on the fact that a 
% multidimensional Gaussian filter can be separated into a series
% of 1D Gaussian filters, applied successively for each dimension.
%
% J=IMGAUSSIAN(I,SIGMA,FILTERSIZE)
%
% inputs,
%   I: The 1D, 2D greyscale/color, or 3D input image with
%           data type Single or Double
%   SIGMA: The sigma of the Gaussian kernel
%   FILTERSIZE: Filter size (single value) (default: sigma*6)
%
% outputs,
%   J: The gaussian filtered image
%
% example,
%   I = im2double(imread('peppers.png'));
%   figure, imshow(imgaussian(I,10));
%
% Function is written by D.Kroon University of Twente (September 2009)
%
%   See also IMFILTER

%% Start function
if(~exist('siz','var')), filterSize=sigma*6; end

if(sigma>0)
    % Make 1D Gaussian kernel
    x=-ceil(filterSize/2):ceil(filterSize/2);
    H = exp(-(x.^2/(2*sigma^2)));
    H = H/sum(H(:));
    
    % Filter each dimension with the 1D Gaussian kernels\
    switch ndims(I)
        case 1
            %... 1D array
            I=imfilter(I,H, 'same' ,'replicate');
        case 2
            %... 2D array
            Hx=reshape(H,[length(H) 1]);
            Hy=reshape(H,[1 length(H)]);
            I=imfilter(imfilter(I,Hx, 'same' ,'replicate'), ...
                Hy, 'same' ,'replicate');
        case 3
            %... 3D array
            if(size(I,3)<4) % Detect if 3D or color image
                Hx=reshape(H,[length(H) 1]);
                Hy=reshape(H,[1 length(H)]);
                for k=1:size(I,3)
                    I(:,:,k)=imfilter(imfilter(I(:,:,k), ...
                        Hx, 'same' ,'replicate'),Hy, 'same' ,'replicate');
                end
            else
                Hx=reshape(H,[length(H) 1 1]);
                Hy=reshape(H,[1 length(H) 1]);
                Hz=reshape(H,[1 1 length(H)]);
                I=imfilter(imfilter(imfilter(I, ...
                    Hx, 'same' ,'replicate'),Hy, 'same' ,'replicate'), ...
                    Hz, 'same' ,'replicate');
            end
        otherwise
            error('imgaussian:input','unsupported input dimension');
    end
end