import pandas as pd
import numpy as np

from skimage import measure
from skimage.filters import threshold_otsu, rank
from skimage.io import imread
from skimage.io import imsave
from skimage import filters
from skimage import morphology
from skimage.data import cells3d
from skimage.measure import label

from pyclesperanto_prototype import imshow
import pyclesperanto_prototype as cle

import gc
from itertools import chain
# define nones functions
nones = lambda n: [None for _ in range(n)]

from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt
import stackview

import sys
#from sys import argv, exit
import os

import pandas as pd
from skimage import measure
from skimage.morphology import disk
from skimage.filters import threshold_otsu, rank
from skimage import measure
from skimage.morphology import disk
from skimage.filters import threshold_otsu, rank
import pyclesperanto_prototype as cle
from skimage import filters
from skimage.filters import try_all_threshold
import napari_simpleitk_image_processing as nsitk


## test https://github.com/haesleinhuepf/spheri/blob/2427867b86b12e814c49ee80a0dc68ac3371b99f/demo/demo.ipynb
## original code doesn't work for some reasons, dissect the functions
DEFAULT_SMOOTH_ITERATIONS = 15

def sphericity_wadell(volume, surface_area):
    # "Hakon Wadell defined sphericity as the surface area of a sphere of the same volume as the particle divided by the actual surface area of the particle."
    # https://en.wikipedia.org/wiki/Sphericity
    # https://www.journals.uchicago.edu/doi/10.1086/624298
    import math
    return pow(math.pi, 1/3) * pow(6 * volume, 2/3) / surface_area

def sphericity_legland(volume, surface_area):
    # https://imagej.net/plugins/morpholibj#shape-factors
    import math
    return 36 * math.pi * pow(volume, 2) / pow(surface_area, 3)

def solidity(surface_volume, convex_hull_volume):
    return surface_volume / convex_hull_volume

def surface_meshes(label_image, smooth_iterations:int=DEFAULT_SMOOTH_ITERATIONS):
    import numpy as np
    import vedo

    result = {}
    for label in np.unique(label_image):
        if label == 0: # skip background
            continue

        binary_image = np.asarray((label_image == label) * 1)

        extended_binary_image = np.zeros([s + 2 for s in binary_image.shape])
        extended_binary_image[1:-1,1:-1,1:-1] = binary_image

        volume = vedo.Volume(extended_binary_image)
        surface = volume.isosurface(value=0.5)

        surface = surface.clean()

        # might work but takes very long time:
        # surface = surface.smooth_mls_2d()
        # might work but takes long time:
        if smooth_iterations > 0:
            surface = surface.smooth(niter=smooth_iterations)

        result[label] = surface

    return result

def soli_measure(label_image, smooth_iterations:int=DEFAULT_SMOOTH_ITERATIONS):
    import numpy as np
    import pandas as pd
    import math
    import vedo

    result = {
        "label":[],
        "surface_area":[],
        "volume":[],
        "convex_hull_area":[],
        "convex_hull_volume":[],
        "solidity":[],
        "sphericity_wadell":[],
        "sphericity_legland":[],
    }

    for label, surface in surface_meshes(label_image, smooth_iterations=smooth_iterations).items():

        convex_hull = vedo.shapes.ConvexHull(surface)

        surface_area = surface.area()
        surface_volume = surface.volume()

        convex_hull_area = convex_hull.area()
        convex_hull_volume = convex_hull.volume()

        result["label"].append(label)
        result["surface_area"].append(surface_area)
        result["volume"].append(surface_volume)
        result["convex_hull_area"].append(convex_hull_area)
        result["convex_hull_volume"].append(convex_hull_volume)
        result["solidity"].append(solidity(surface_volume, convex_hull_volume))
        result["sphericity_wadell"].append(sphericity_wadell(surface_volume, surface_area))
        result["sphericity_legland"].append(sphericity_legland(surface_volume, surface_area))

    return pd.DataFrame(result)



if __name__ == "__main__":
    image = sys.argv[1]
    #num1 = sys.argv[1]
    #num2 = sys.argv[2]
    #num3 = sys.argv[3]

    imageDir =os.path.dirname(image)
    imageName = os.path.basename(image)
    fileName = imageName.replace('_blur_cp_masks_cleaned.tif','')

    print(fileName)

    outDir = '/groups/tanaka/People/current/jiwang/projects/RA_competence/images_data/results/' + 'embryo_features' + '/'
    if not os.path.exists(outDir):
        os.mkdir(outDir)
        #print(outDir)

    mask = imread(os.path.join(imageDir, str(fileName + "_blur_cp_masks_cleaned.tif"))) # mask segmeted cysts
    crops = imread(os.path.join(imageDir, str(fileName + ".tif")))

    foxa2 = crops[:, :, :, 0]
    pax6 = crops[:, :, :, 1]
    sox2 = crops[:, :, :, 2]
    dapi = crops[:, :, :, 3]

    labels_mask, nb_mask = measure.label(mask, return_num = True)

    print(mask.shape, crops.shape, nb_mask)

    ## sphericity parameters
    sphere_stats = soli_measure(labels_mask)


    ## mask and foxa2 channel features
    props_C1 = measure.regionprops_table(labels_mask,
            intensity_image = foxa2,
            properties=('label', 'area', 'equivalent_diameter_area', 'centroid',
            'intensity_max', 'intensity_mean', 'intensity_min', 'intensity_std'))
    props_C1 = pd.DataFrame(props_C1)
    props_C1 = props_C1.add_suffix('_foxa2')

    ## C2
    props_C2 = measure.regionprops_table(labels_mask,
                        intensity_image = pax6,
                        properties=('intensity_max', 'intensity_mean', 'intensity_min', 'intensity_std'))
    props_C2 = pd.DataFrame(props_C2)
    props_C2 = props_C2.add_suffix('_pax6')

    # channel C3
    props_C3 = measure.regionprops_table(
            labels_mask,
            intensity_image = sox2,
            properties=('intensity_max', 'intensity_mean', 'intensity_min', 'intensity_std'))
    props_C3 = pd.DataFrame(props_C3)
    props_C3 = props_C3.add_suffix('_sox2')

    # chanel 4
    props_C4 = measure.regionprops_table(
            labels_mask,
            intensity_image = dapi,
            properties=('intensity_max', 'intensity_mean', 'intensity_min', 'intensity_std'))
    props_C4 = pd.DataFrame(props_C4)
    props_C4 = props_C4.add_suffix('_dapi')


    #df = pd.concat([props_C1, props_C2, props_C3, props_C4], axis=1)
    df = pd.concat([sphere_stats, props_C1, props_C2, props_C3, props_C4], axis=1)
    
    df.to_csv(os.path.join(outDir, "featuresCollection_" + fileName + ".csv"), index=True, header=True)
