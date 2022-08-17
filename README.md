# fishSearch
## Purpose 
This module makes usage of various image segmenation algorithms to analysis smFISH data. This guide will detail steps after acquistion pertaining to probe analysis. 

## Installation 
First on your local machine use git clone to create your version of the repository. 

### Example: git clone
```
git clone https://github.com/choulab/fishSearch.git
```

### Example: update clone
```
#cd into working directory
git pull
```

## Create Conda Environment
Ensure that you have conda installed onto your computer and is included in path. Then open terminal/command line is the newly downloaded repository. In the repository use the 'environment.yml' file to create a new conda environment.

### Example: create Conda environment
```
conda env create --name envname --file=environment.yml
```

### Example: heck environment was created
```
conda info --envs
```

## Example Workflow:

### Open Jupyter Lab in Cloned Directory
```
jupyter lab
```

### Start Analysis
#### Import all Necessary Functions
``` python
from fishSearch.fishSearch import *
```

#### Read in Datafile
``` python
#function reads in img file and returns datastruct and import file anme
img, img_name = read_img()
```

#### Seperate Channels
``` python
#function takes in tiff object and returns seperated images based on number of channels
seperate_channels(img, img_name)
```

#### Generate Views of Channel
``` python
#function takes in img object and creates 2D and 3D projection of the image stack
views = generate_projects(img, img_name)
DAPI_3D = views[0]
DAPI_2D = views[1]
RNA_3D = views[2]
RNA_2D = views[3]
```

#### Create nucleus and cell segmentation masks
``` python
#function takes in max projection of DAPI channel and both RNA channels - diameter and threshold are last two variables respectively
nuc_label, cell_label = cell_nuc_segmentation(DAPI_2D, RNA_3D, RNA_2D, 100, 40)
```

#### Take in RSFISH Spot Detection and Map Detection on Max RNA Projection
``` python
#function takes in Max RNA projection object and takes input for RNA detection for RSFISH. Returns object for mapped RNA
rs_spots = rsfish_analysis(RNA_2D)
```

#### Extract Cell Specific Information 
``` python
#function takes in nucleus mask, cell mask, mapped RNA spots, and Max RNA Projection. Returns object that contains cell specific information. 
fov_results = cell_extraction(nuc_label, cell_label, rs_spots, RNA_2D)
```

#### Visualize Detected Cells and Internal Spots
``` python
#function takes in object with cell specific information and provides all detected cells and mapped spots
cell_level_visualization(fov_results)
```

#### Create CSV in Data Folder under Graph Data for Average Spot Count, maxMinCount, and Spot Intensity per Cell
``` python
spotCount(fov_results,img_name)
maxMinCount(fov_results,img_name)
spotIntensity(cell_label, rs_spots, img_name)
```
