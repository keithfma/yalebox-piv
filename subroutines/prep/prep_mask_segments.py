
import argparse
from skimage.segmentation import felzenszwalb
from scipy.io import loadmat, savemat

ap = argparse.ArgumentParser(
    description='Multiband Felzenswalb image segmentation',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
ap.add_argument('input_file', type=str,
                help='Input MAT file containing multiband image array')
ap.add_argument('output_file', type=str,
                help='Output MAT file to save segment labels (superpixels)')
ap.add_argument('--scale', type=int, required=False, default=100,
                help='TODO')
ap.add_argument('--sigma', type=float, required=False, default=0.5,
                help='TODO')
ap.add_argument('--min_size', type=int, required=False, default=50,
                help='TODO')
args = ap.parse_args()

input_dict = loadmat(file_name=args.input_file)
input_var_names = [k for k in input_dict.keys() if not k.startswith('__')]
assert len(input_var_names) == 1, 'Too many variables in input mat file'
input_data = input_dict[input_var_names[0]]

segments = felzenszwalb(
    input_data, scale=args.scale, sigma=args.sigma, min_size=args.min_size, multichannel=True)

output_dict = {'segments': segments}
savemat(file_name=args.output_file, mdict=output_dict)


