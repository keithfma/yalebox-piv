function [] = movie_single(input_file, output_file, varargin)
%
% Required Arguments:
%   input_file: String, path to input netCDF file as produced by prep_series()
%   output_file: String, path to output video file, without file extension
%
% Optional Parameters ('Name', Value)::
%   'VideoProfile': String, video file type, see MATLAB VideoWriter help for
%       options, default = 'MPEG-4'
%   'FrameRate': Integer, video frames-per-second, default = 10
%   'Color': Logical, set True to use RGB image series for movie, or False to
%       use (equalized) grayscale image series, default = True
%   'Mask': Logical: set True to apply mask to images, or False to use the
%       images as-is, default = True
%   'MinThresh': Minumum intensity value, pixels dimmer than this are set to 0,
%       must be in the range [0, 1], default = 0
%   'MaxDim': Vector, [width, height] maximum dimensions of the output video,
%       image (and coordinates) are resized to fit in this bounding box, note
%       that MATLAB places some limits on this depending on your screen size,
%       default = [1080, 1920]
%   'SptPosition': Tip location of s-point triangles in [x, y] world coordinates as
%       (number of triangles) x 2 matrix, with 1 row for each triangle,
%       default = []
%   'SptLength': Scalar, side length for all (equilateral) s-point triangles in
%       world coordinates, default = 0.01
%   'SptColor': String, s-point triangel color, as a color definition MATLABcan
%       understand, default = 'r'
%   'SptAlpha': Scalar, opacity (a.k.a. alpha) for the s-point triangles, must
%       be in the range [0, 1], default = 1
%   'TitleStr': String, or cell array of strings for each frame in the movie,
%       default = '' (no title)
%   'TitleSize': Scalar, integer, title font size in points, default = 18
%   'TitleColor': Title color, as any definition MATLAB can understand,
%       default = 'r'
%   'TitleBoxColor: Title background box color, as any definition MATLAB can
%       understand, default = 'w'
%   'TitleBoxAlpha': Scalar, opacity (a.k.a. alpha) for the title background
%       box, must be in the range [0, 1], default = 0
%   'TitleLocation': String, title horizontal location, must be one of
%       {'left', 'center', 'right'}, default = 'center'
%   'CounterPosition': Counter position in [x, y] world coordinates,
%       default = [] (no counter)
%   'CounterSize': Scalar, integer, counter font size in points, default = 18
%   'CounterColor: Counter text color, as any definition MATLAB can
%       understand, default = 'r'
%   'CounterBoxColor: Counter background box color, as any definition MATLAB can
%       understand, default = 'w'
%   'CounterBoxAlpha': Scalar, opacity (a.k.a. alpha) for the counter background
%       box, must be in the range [0, 1], default = 0 (no box)
%   'ScalePosition: Vector, scalebar location and size as a MATLAB position
%       vector in world units, [left, bottom, width, height], default = []
%       (no scale bar)
%   'ScaleColor': Scalebar fill color, as any definition MATLAB can understand,
%       default = 'r'
%   'ScaleAlpha': Scalar, opacity (a.k.a. alpha) for the scalebar fill color
%       must be in the range [0, 1], default = 1 (solid)
%   'ScaleLabel': String, scalebar label text, default = ''
%   'ScaleLabelSize': Scalar, integer, scalebar label font size in points,
%       default = 18 
%   'ScaleLabelColor': Scalebar label color, as any definition MATLAB can
%       understand, default = 'r'
%   'ScaleBoxColor': Scalebar background box color, as any definition MATLAB can
%       understand, default = 'w'
%   'ScaleBoxAlpha': Scalar, opacity (a.k.a. alpha) for the counter background
%       box, must be in the range [0, 1], default = 0 (no box)
%   'TestFrame': Integer, display the the frame with the given index and do not
%       make any video, intended for parameter exploration, default = 0 (off)



%% parse inputs

parser = inputParser();

parser.addRequired('input_file', @(x) exist(x, 'file') == 2);
parser.addRequired('output_file', @ischar);
parser.addParameter('VideoProfile', 'MPEG-4', @ischar);
parser.addParameter('FrameRate', 10, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer'}) ); 
parser.addParameter('Color', true, ...
    @(x) validateattributes(x, {'numeric', 'logical'}, {'scalar', 'binary'}) ); 
parser.addParameter('Mask', true, ...
    @(x) validateattributes(x, {'numeric', 'logical'}, {'scalar', 'binary'}) ); 
parser.addParameter('MinThresh', 0, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', '>=', 0, '<=', 1}) ); 
parser.addParameter('MaxDim', [1920, 1080], ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', 2, 'positive', 'integer'}) );
parser.addParameter('SptPosition', [], ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}) );
parser.addParameter('SptLength', 0.01, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}) );
parser.addParameter('SptColor', 'r', ...
    @(x) validateattributes(x, {'numeric', 'char'}, {}) );
parser.addParameter('SptAlpha', 1, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', '>=', 0, '<=', 1}) );
parser.addParameter('TitleStr', '', ...
    @(x) validateattributes(x, {'char', 'cell'}, {}) );
parser.addParameter('TitleSize', 18, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', 'positive'}) );
parser.addParameter('TitleColor', 'r', ...
    @(x) validateattributes(x, {'numeric', 'char'}, {}) );
parser.addParameter('TitleBoxColor', 'w', ...
    @(x) validateattributes(x, {'numeric', 'char'}, {}) );
parser.addParameter('TitleBoxAlpha', 0, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', '>=', 0, '<=', 1}) );
parser.addParameter('TitleLocation', 'center', ...
    @(x) validateattributes(x, {'char'}, {'nonempty'}) );
parser.addParameter('CounterPosition', [], ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}) );
parser.addParameter('CounterSize', 18, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', 'positive'}) );
parser.addParameter('CounterColor', 'r', ...
    @(x) validateattributes(x, {'numeric', 'char'}, {}) );
parser.addParameter('CounterBoxColor', 'w', ...
    @(x) validateattributes(x, {'numeric', 'char'}, {}) );
parser.addParameter('CounterBoxAlpha', 0, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', '>=', 0, '<=', 1}) );
parser.addParameter('ScalePosition', [], ...
    @(x) validateattributes(x, {'numeric'}, {'1d'}) );
parser.addParameter('ScaleColor', 'r', ...
    @(x) validateattributes(x, {'numeric', 'char'}, {}) );
parser.addParameter('ScaleAlpha', 1, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', '>=', 0, '<=', 1}) );
parser.addParameter('ScaleLabel', '', ...
    @(x) validateattributes(x, {'char'}, {}) );
parser.addParameter('ScaleLabelSize', 18, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', 'positive'}) );
parser.addParameter('ScaleLabelColor', 'r', ...
    @(x) validateattributes(x, {'numeric', 'char'}, {}) );
parser.addParameter('ScaleBoxColor', 'w', ...
    @(x) validateattributes(x, {'numeric', 'char'}, {}) );
parser.addParameter('ScaleBoxAlpha', 1, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', '>=', 0, '<=', 1}) );
parser.addParameter('TestFrame', 0, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer'}) ); 

parser.parse(input_file, output_file, varargin{:});
args = parser.Results;

% define mode
movie_mode = ~args.TestFrame;

% warn about unusued parameters
if args.MinThresh ~= 0 && args.Color
    warning('MinThresh is ignored for Color images');
end

% expand title strings to cell array
ns = length(ncread(args.input_file, 'step'));
if ischar(args.TitleStr)
    constant_title = args.TitleStr;
    args.TitleStr = cell(ns, 1);
    args.TitleStr(:) = {constant_title};
elseif iscell(args.TitleStr)
    assert(length(args.TitleStr == ns));
end
    
%% Create video

% get coordinate vars
xx = double(ncread(args.input_file, 'x'));
yy = double(ncread(args.input_file, 'y'));
step = double(ncread(args.input_file, 'step'));

% prepare movie object
if movie_mode
    writer = VideoWriter(args.output_file, args.VideoProfile);
    writer.FrameRate = args.FrameRate;
    writer.open();
end

% make frames
for ii = 1:numel(step)
    
    if ~movie_mode && step(ii) ~= args.TestFrame
        % skip all but indicated frame
        continue
    end
    
    fprintf('step: %i\n', step(ii));
    
    % prepare image
    if args.Color
        frame = ncread(input_file, 'img_rgb', [1, 1, 1, ii], [inf, inf, inf, 1]);
    else
        frame = ncread(input_file, 'img', [1, 1, ii], [inf, inf, 1]);
    end
    frame = im2double(frame);
    
    if args.Mask
        mask_manual = ncread(input_file, 'mask_manual');
        mask_auto = ncread(input_file, 'mask_auto', [1, 1, ii], [inf, inf, 1]);
        mask = mask_manual & mask_auto;
        frame(repmat(~mask, [1, 1, size(frame, 3)])) = 0;
    end
    
    [frame, xf, yf] = movie_frame_resize(frame, xx, yy, args.MaxDim);
    
    [frame, yf] = movie_frame_flip(frame, yf);
    
    if args.MinThresh > 0 && ~args.Color
        frame(frame < args.MinThresh) = 0;
    end
    
    fprintf('%d x %d\n', size(frame, 2), size(frame, 1));
    
    % annotate image
    frame = movie_frame_spoint(frame, xf, yf, args.SptPosition, ...
        args.SptLength, args.SptColor, args.SptAlpha);
    
    frame = movie_frame_title(frame, args.TitleStr{ii}, args.TitleSize, ...
        args.TitleColor, args.TitleBoxColor, args.TitleBoxAlpha, ...
        args.TitleLocation);
    
    if ~isempty(args.CounterPosition)
        frame = movie_frame_counter(frame, xf, yf, step(ii), ...
            args.CounterPosition, args.CounterSize, args.CounterColor, ...
            args.CounterBoxColor, args.CounterBoxAlpha);
    end
    
    if ~isempty(args.ScalePosition)
        frame = movie_frame_scalebar(frame, xf, yf, args.ScaleLabel, ...
            args.ScalePosition, args.ScaleColor, args.ScaleAlpha, ...
            args.ScaleLabelSize, args.ScaleLabelColor, args.ScaleBoxColor, ...
            args.ScaleBoxAlpha);
    end
    
    % add to movie writer
    if movie_mode
        writer.writeVideo(frame);
    end

end

if movie_mode
    writer.close();
else
    hf = figure;
    hf.Name = sprintf('Movie Frame for Step %i', args.TestFrame);
    imshow(frame);
end