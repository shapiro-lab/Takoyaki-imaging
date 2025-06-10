import os
import scipy.io
import numpy as np
import skimage.io
import napari
#from naparimovie import Movie
from napari_animation import Animation
#gui qt5

script_dir = os.path.dirname(os.path.abspath(__file__))

def update_slider(event):

    time = viewer.dims.current_step[0]
    if time < 11:
       viewer.text_overlay.text = ''
    if time >= 11 and time < 28:
       viewer.text_overlay.text = 'Being injected'
    if time >= 28:
       viewer.text_overlay.text = ''

subfolder_name = 'TakoyakiAM_live_imaging'
frame_num = 150

mat = np.zeros(shape = (frame_num, 171, 153, 136)) # 171, 153, 133

for i in range(1, frame_num+1):
    data_file_name = f'pAM3D_{i}.mat'

    data_file_path = os.path.join(script_dir, subfolder_name, data_file_name)

    import_file = scipy.io.loadmat(data_file_path)
    mat[i-1,] = import_file[f'pAM3D_{i}'][:, ::-1, :]


Dop_file_path = os.path.join(script_dir, f'Dop.mat')
Dop_mat = scipy.io.loadmat(Dop_file_path)[f'Dop_final'] # size 98 98 102

viewer = napari.view_image(mat[:, 85:, :, :128], colormap = 'green', gamma = 1.0, contrast_limits = [15000, 100000])
viewer.dims.ndisplay = 3
viewer.camera.angles = (-90.0, 0.0, 0.0) # (-91.89697387647686, -41.77534611484706, 138.3101242566956)
viewer.camera.zoom = 2
viewer.axes.visible = True
viewer.text_overlay.visible = True
viewer.dims.events.current_step.connect(update_slider)

viewer.add_image(Dop_mat[49:,::-1,33:98], colormap = 'red', gamma = 1.10, contrast_limits = [160000, 1456962], 
                 blending = 'additive', scale = [2, 2, 203/102], translate = [0.5, -21, 0]) # -12 -21 -66


# Add a white border around the volume
cropped = Dop_mat[49:,::-1,33:98]
corners = np.array([
    [0, 0, 0],
    [0, cropped.shape[1]*2, 0],
    [cropped.shape[0]*2, cropped.shape[1]*2, 0], 
    [cropped.shape[0]*2, 0, 0],
    [0, 0, cropped.shape[2]*203/102],
    [0, cropped.shape[1]*2, cropped.shape[2]*203/102],
    [cropped.shape[0]*2, cropped.shape[1]*2, cropped.shape[2]*203/102],
    [cropped.shape[0]*2, 0, cropped.shape[2]*203/102]
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
    translate = [0, -21, 0]
)

napari.run()
