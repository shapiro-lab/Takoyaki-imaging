import os
import scipy.io
import numpy as np
import skimage.io
import napari
#from naparimovie import Movie
from napari_animation import Animation
#gui qt5

script_dir = os.path.dirname(os.path.abspath(__file__))

# Add pAM image in napari
frame_num = 1

mat = np.zeros(shape = (frame_num, 171, 153, 136)) # 171, 153, 136

for i in range(1, frame_num+1):
    data_file_name = f'TakoyakiAM.mat'

    data_file_path = os.path.join(script_dir, data_file_name)

    import_file = scipy.io.loadmat(data_file_path)
    mat = import_file[f'pAM3D_{i}'][:, ::-1, :]

viewer = napari.view_image(mat[:86, :, 2:122], colormap = 'green', gamma = 2.0, translate = (12, 21, 0), contrast_limits = [10000, 31127])
viewer.dims.ndisplay = 3
viewer.camera.angles = (-90.0, 0.0, 90.0) # [-90, -27, 125]
viewer.camera.zoom = 2
viewer.axes.visible = True

# Add doppler image in napari
Dop_file_path = os.path.join(script_dir, f'Dop.mat')
Dop_mat = scipy.io.loadmat(Dop_file_path)[f'Dop_final'] # size 98 98 102
Dop_mat = Dop_mat[:, ::-1, :]

viewer.add_image(Dop_mat[6:49,11:-11,34:94], colormap = 'red', gamma = 1.30, contrast_limits = [36450, 712340], 
                 blending = 'additive', scale = [2, 2, 2], translate = [12, 22, 0.5]) # -12 -21 -66 [33:99] 203/102

# Add BURST image in napari
BURST_file_path = os.path.join(script_dir, f'BURST_braintumor.mat')
BURST_mat = scipy.io.loadmat(BURST_file_path)[f'svd_filtered_final'] # size 98 98 102

viewer.add_image(BURST_mat[:86,::-1,2:122], colormap = 'cyan', gamma = 1, contrast_limits = [10000, 120000], 
                 blending = 'additive', translate = [12, 21, 0]) # [30, 42, 28]

# Add a white border around the volume
cropped = mat[:86, :, 2:122]
corners = np.array([
    [0, 0, 0],
    [0, cropped.shape[1], 0],
    [cropped.shape[0], cropped.shape[1], 0],
    [cropped.shape[0], 0, 0],
    [0, 0, cropped.shape[2]],
    [0, cropped.shape[1], cropped.shape[2]],
    [cropped.shape[0], cropped.shape[1], cropped.shape[2]],
    [cropped.shape[0], 0, cropped.shape[2]]
])

lines = np.array([
    [corners[0], corners[1]],
    [corners[1], corners[2]],
    [corners[2], corners[3]],
    [corners[3], corners[0]],
    [corners[4], corners[5]],
    [corners[5], corners[6]],
    [corners[6], corners[7]],
    [corners[7], corners[4]],
    [corners[0], corners[4]],
    [corners[1], corners[5]],
    [corners[2], corners[6]],
    [corners[3], corners[7]]
])

viewer.add_shapes(
    lines,
    shape_type='line',
    edge_color='white',
    edge_width=0.5,
    name='Volume Border',
    translate= [12, 21, 0]
)

napari.run()

