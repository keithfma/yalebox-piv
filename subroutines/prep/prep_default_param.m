function pp = prep_default_param()
% function pp = prep_default_param()
%
% Return parameter struct for image pre-processing with some default values
% %

% test ----------

pp.test.num_images.help = "int, number of images to include for test run";
pp.test.num_images.value = 20;

pp.test.test_file.help = "string, step image filename (without directory) to use for testing"; 
pp.test.test_file.value = "";

% images ---------

pp.images.path.help = "string, directory containing all image files";
pp.images.path.value = ""; 

pp.images.woco_file.help = "string, world coordinate image filename (without directory)";
pp.images.woco_file.value = "";

pp.images.exp_files.help = "list of strings, all experiment image filenames (without directory), assumed to be ordered by time";
pp.images.exp_files.value = [];

% woco ----------

% 		
% 	"woco": {
% 		"xw": {
% 			"help": "1D vector, x-position of control points in world coordinates",
% 			"value": [-0.2,-0.2,-0.2,-0.1,-0.1,-0.1,0,0,0,0.1,0.1,0.1,0.2,0.2,0.2,0.3,0.3,0.3,0.4,0.4,0.4,0.5,0.5,0.5,0.6,0.6]
% 		},
% 		"yw": {
% 			"help": "1D vector, x-position of control points in world coordinates",
% 			"value": [0,0.1,0.2,0,0.1,0.2,0,0.1,0.2,0,0.1,0.2,0,0.1,0.2,0,0.1,0.2,0,0.1,0.2,0,0.1,0.2,0,0.1]
% 		},
% 		"xp": {
% 			"help": "1D vector, x-position of control points in pixel coordinates",
% 			"value": [147.6801356,147.6265596,148.1573209,660.2430281,659.2737798,658.8500566,1171.390491,1170.568882,1170.850133,1689.332289,1687.850734,1687.609345,2206.411138,2205.67073,2204.629877,2723.928475,2722.370148,2720.591839,3238.277484,3237.009267,3235.030908,3749.723177,3747.047225,3743.75119,4260.407055,4257.238727]
% 		},
% 		"yp": {
% 			"help": "1D vector, x-position of control points in pixel coordinates",
% 			"value": [1398.529661,874.8481166,362.8539454,1393.489722,869.9479869,358.5694692,1390.863229,865.8081966,354.5965107,1386.573585,861.5229822,350.9800998,1385.550449,863.6349544,349.1316628,1384.536909,861.6425696,348.7655187,1380.572809,861.7850409,350.0765507,1382.480301,860.8091489,351.2975231,1381.155244,861.3148116]
% 		}
% 	},
% 	"crop": {
% 		"xlim": {
% 			"help": "2-element array, minimum and maximum x coordinates for crop in world coordinate units (m)",
% 			"value": [-0.25,0.6]
% 		},
% 		"ylim": {
% 			"help": "2-element array, minimum and maximum y coordinates for crop in world coordinate units (m)",
% 			"value": [0.001,0.2]
% 		},
% 		"npts": {
% 			"help": "int, number of points to use in local weighted mean calculation for image warping",
% 			"value": 10
% 		}
% 	},
% 	"mask_manual": {
% 		"poly": {
% 			"help": "2D array, vertices of mask polygons, x-coords in row 1 and y-coords in row 2, polygons separated by NaN",
% 			"value": [
% 				[3322.455461,3434.795699,3695.220795,3996.496886,4274.794293,4277.34748,3230.540722,"_NaN_",1032.216363,1093.763148,1228.511757,1030.171619,1030.171619,"_NaN_"],
% 				[688.0635434,499.1276895,361.2555799,328.0641461,325.5109589,1048.062941,1040.403379,"_NaN_",44.04107306,47.10818859,0.07908374131,0.07908374131,44.45002179,"_NaN_"]
% 			]
% 		}
% 	},
% 	"mask_auto": {
% 		"hue_lim": {
% 			"help": "2-element vector, double, range [0, 1]. [minimum, maximum] HSV 'hue' included as sand in the mask",
% 			"value": [0.01,0.3]
% 		},
% 		"value_lim": {
% 			"help": "2-element vector, double, range [0,1]. [minimum, maximum] HSV 'value' included as sand in the mask",
% 			"value": [0.07,0.6]
% 		},
% 		"entropy_lim": {
% 			"help": "2-element vector, double, range [0, 1]. [minimum, maximum] entropy included as sand in the mask",
% 			"value": [0.55,1]
% 		},
% 		"entropy_len": {
% 			"help": "scalar, integer, window size in pixels for entropy filter",
% 			"value": 11
% 		},
% 		"morph_open_rad": {
% 			"help": "scalar, double, radius of disk structuring element used in mophological opening filter",
% 			"value": 10
% 		},
% 		"morph_erode_rad": {
% 			"help": "scalar, double, radius of disk structuring element used in mophological erosion filter",
% 			"value": 3
% 		}
% 	},
% 	"intensity": {
% 		"eql_len": {
% 			"help": "scalar, integer, odd. Side length (in pixels) for the local neighborhood used to compute the transform for each pixel",
% 			"value": 41
% 		}
% 	}
% }
