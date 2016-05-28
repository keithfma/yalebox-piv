% Parameters for single-view movie of fault_ss_02_siden - RGB & streak  
% %

%% Define parameters

% source parameters
prm.input_file = '~/Documents/dissertation/yalebox-exp-fault/data/fault_ss_02_siden_image.nc'; 
prm.output_stub = '~/Documents/dissertation/yalebox-exp-fault/movie/fault_ss_02_siden_movie_single'; 

% video parameters
prm.frame_rate = 8;
prm.max_dim = [1920, 1080];
prm.streak_threshold = 0.95;  
prm.streak_memory = 0.5; 

% s-point triangle annotation parameters
prm.tri_tip = [0, 0; -0.15, 0]; % 1 triangle per row, position in m
prm.tri_len = 0.01; % m
prm.tri_color = 'red';
prm.tri_opacity = 1;

% title annotation parameters
prm.title_str = {'Loose Sand', 'Sieved Sand'};
prm.title_str_start = [0, 227]; 
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

show_frame_number = 227;
movie_single(prm, 'color', show_frame_number);
movie_single(prm, 'streak', show_frame_number);

%% Save final parameters

param_name = [prm.output_stub, '_param.mat'];
fprintf('saving parameters as: %s\n', param_name);
save(param_name, 'prm');

