import numpy as np
import pandas as pd
import seaborn as sns
import nibabel as nib
from scipy.spatial.distance import cdist
from scipy.ndimage.morphology import binary_erosion, binary_dilation

def find_border(data):
    eroded = binary_erosion(data)
    border = np.logical_and(data, np.logical_not(eroded))
    return border

def get_coordinates(data, affine):
    if len(data.shape) == 4:
        data = data[:, :, :, 0]
    indices = np.vstack(np.nonzero(data))
    indices = np.vstack((indices, np.ones(indices.shape[1])))
    coordinates = np.dot(affine, indices)
    return coordinates[:3, :]

def save_nifti(img,affine,header,filename):
    out = nib.Nifti1Image(img, affine, header)
    nib.save(out,filename)

# Load images
gm                   = nib.load(snakemake.input.ribbon[0])
vessels              = nib.load(snakemake.input.vessels)
diameters            = nib.load(snakemake.input.diameter)

# Identify vessel borders
vesselsBorder        = find_border(vessels.get_fdata())
vesselsBorder_coords = get_coordinates(vesselsBorder, vessels.affine)

# Identify GM tissue
gm_data              = gm.get_fdata()
gm_data[gm_data!=1]  = 0
gm_data              = binary_dilation(gm_data, iterations=10)
gm_coords            = get_coordinates(gm_data, gm.affine)

# Calculate distance from each GM voxel to closest vessel and its diameter
min_distance         = np.zeros((gm_coords.shape[1],1))
closest_diameter     = np.zeros((gm_coords.shape[1]))

for i in range(0,gm_coords.shape[1]):
    dist                = None
    dist                = cdist(gm_coords.T[i].reshape(1,3), vesselsBorder_coords.T)
    min_distance[i,0]   = np.amin(dist)
    closest_diameter[i] = np.argmin(dist)

vessels = vessels.get_fdata()
vessels = vessels.flatten()    

img             = np.zeros(gm.get_fdata().shape)
closestDiameter = np.zeros(gm.get_fdata().shape)

diameters = diameters.get_fdata()
diameters = diameters[np.nonzero(vesselsBorder)]
diameters = diameters.flatten()

# Populate new nifti images (only within GM voxels)
indices = np.nonzero(gm_data)
for i in range(0,len(min_distance)):
    x = indices[0][i]
    y = indices[1][i]
    z = indices[2][i]
    img[x,y,z]             = min_distance[i]
    closestDiameter[x,y,z] = diameters[int(closest_diameter[i])]

# Save nifti images
save_nifti(img, gm.affine, gm.header, snakemake.output.distance)
save_nifti(closestDiameter, gm.affine, gm.header, snakemake.output.diameter)
save_nifti(vesselsBorder, gm.affine, gm.header, snakemake.output.edges)