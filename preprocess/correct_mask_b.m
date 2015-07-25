% New version of masking routine, very rough

% read in data

% convert to hsv and separate out the hue and value layers

% convert to grayscale and run through an entroy filter

% run everything through median filters

% create mask by specifying ranges for each: 
%   low hue
%   low value
%   high entropy

% add 1's at sides and bottom
% fill holes
% remove 1's at sides and bottom

% morphological opening

% done!