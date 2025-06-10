import os
import scipy.io
import numpy as np
import skimage.io
import napari
#from naparimovie import Movie
from napari_animation import Animation
#gui qt5

script_dir = os.path.dirname(os.path.abspath(__file__))

frame_num = 1

mat = np.zeros(shape = (frame_num, 171, 153, 136)) # 171, 153, 133

for i in range(1, frame_num+1):
    data_file_name = f'TakoyakiAM.mat'

    data_file_path = os.path.join(script_dir, data_file_name)

    import_file = scipy.io.loadmat(data_file_path)
    mat[i-1,] = import_file[f'pAM3D_150'][:, ::-1, :]


BURST_file_path = os.path.join(script_dir, f'BURST_ventricle.mat')
BURST_mat = scipy.io.loadmat(BURST_file_path)[f'svd_filtered_final'] # size 98 98 102

viewer = napari.view_image(mat[:, 85:, :, :128], colormap = 'green', gamma = 1.0, contrast_limits = [5000, 100000])
viewer.dims.ndisplay = 3
viewer.camera.angles = (-90.0, 0.0, 0.0) # -90.82878631143305, -12.266562341867669, 129.20178889568152
# (-88.83748074013091, -25.17316862723227, -44.21402370198491)
viewer.camera.zoom = 2
viewer.axes.visible = True
viewer.text_overlay.visible = True

viewer.add_image(BURST_mat[85:, ::-1, :128], colormap = 'cyan', gamma = 1.10, blending = 'additive', contrast_limits = [0, 400000]) # -12 -21 -66

napari.run()
