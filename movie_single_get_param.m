% Template script for exploring and selecting parameters for single-view movies
% (i.e. input arguments for movie_single).The intended wusage is to make a copy
% of the script for a given experiment, run it cell by cell, modifying the
% default parameters to suit the experiment particulars. 
% %

%% Define parameters

% source parameters
prm.input_file = '~/Documents/dissertation/yalebox-exp-fault/data/fault_ss_02_siden_image.nc'; 
prm.output_stub = 'fault_ss_02_siden_movie_bw'; 
prm.input_var = 'img'; 

% video parameters
prm.frame_rate = 7;
prm.max_dim = [1920, 1080];
prm.threshold = -inf;
prm.memory = 0;

% s-point triangle annotation parameters
prm.tri_tip = [0, 0; -0.15, 0]; % 1 triangle per row, position in m
prm.tri_len = 0.01; % m
prm.tri_color = 'red';
prm.tri_opacity = 1;

% title annotation parameters
prm.title_str = {'Loose Sand', 'Compacted Sand', 'Sieved Sand'};
prm.title_str_start = [0, 146, 341]; 
prm.title_size = 72;
prm.title_color = 'red';
prm.title_box_color = 'white';
prm.title_box_opacity = 0;

% scalebar annotation parameters
prm.scale_pos = [-0.15, 0.18, 0.1, 0.008]; % [x, y, width, height], in world coords
prm.scale_label = '10 cm';
prm.scale_color = 'red';
prm.scale_opacity = 1;
prm.scale_text_size = 36;
prm.scale_text_color = 'red';
prm.scale_box_color = 'white';
prm.scale_box_opacity = 0;

% counter annotation parameters
prm.count_pos = [0.5, 0.18];
prm.count_size = 36;
prm.count_color = 'red';
prm.count_box_color = 'white';
prm.count_box_opac = 0;

%% Test parameters

show_frame_number = 1;
movie_single(prm, show_frame_number)

