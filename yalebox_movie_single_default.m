function prm = yalebox_movie_single_default()
%
% Creates a struct containing all the parameters needed by
% yalebox_movie_single(), populated with default values. Typical usage
% would be to call this function to get a properly formated parameter
% struct, edit the parameters as needed, then pass the edited struct to
% yalebox_movie_single() to make the movie.
% 
% Arguments:
%   
%   prm = Struct, see yalebox_movie_single help for more information
%   
% Keith Ma, August 2015

% video
prm.frame_rate = 5;
prm.max_dim = [1920, 1080];
prm.threshold = 0;
prm.decay_factor = 0;

% s-point triangle annotation
prm.tri_tip = [0, 0; -0.15, 0]; % 1 triangle per row, position in m
prm.tri_len = 0.01; % m
prm.tri_color = 'red';
prm.tri_opacity = 0.7;

% title annotation
prm.title_str = {'Loose Sand', 'Compacted Sand', 'Sieved Sand'};
prm.title_str_start = [0, 146, 341]; 
prm.title_size = 72;
prm.title_color = 'red';
prm.title_box_color = 'white';
prm.title_box_opacity = 1;

% scalebar annotation
prm.scale_pos = [-0.1, 0.15, 0.1, 0.01]; % [x, y, width, height], in world coords
prm.scale_label = '10 cm';
prm.scale_color = 'red';
prm.scale_opacity = 1;
prm.scale_text_size = 72;
prm.scale_text_color = 'red';
prm.scale_box_color = 'white';
prm.scale_box_opacity = 1;

% counter annotation
prm.count_pos = [0.3, 0.15];
prm.count_size = 72;
prm.count_color = 'red';
prm.count_box_color = 'white';
prm.count_box_opac = 1;
