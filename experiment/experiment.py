"""
Image segmentation comparison, see:
http://scikit-image.org/docs/dev/auto_examples/segmentation/plot_segmentations.html
"""

import matplotlib.pyplot as plt
import numpy as np

from skimage.util import img_as_float
from skimage.io import imread
from skimage.segmentation import felzenszwalb, mark_boundaries
from skimage.color import rgb2hsv
from skimage.filters.rank import entropy
from skimage.morphology import disk

# constants
img_file = 'fault_ss_01_siden_250.jpg'

# load image
img = img_as_float(imread(img_file))

# transform image
# TODO: change color space, add entropy filters, other?
hsv = rgb2hsv(img)
entropy_kernel = disk(5);
h_e = entropy(hsv[:, :, 0], entropy_kernel)
s_e = entropy(hsv[:, :, 1], entropy_kernel)
v_e = entropy(hsv[:, :, 2], entropy_kernel)

img_plus = np.dstack((hsv, h_e, s_e, v_e))

# display
plt.imshow(v_e)
plt.tight_layout()
plt.show()

# segment
segments = felzenszwalb(img_plus, scale=500, sigma=0.5, min_size=50, multichannel=True)

# display
plt.imshow(mark_boundaries(img, segments))
plt.tight_layout()
plt.show()

