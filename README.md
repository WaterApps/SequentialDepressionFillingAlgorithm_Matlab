# Sequential Depression Filling Algorithm

## Getting started:
Currently we support geotiff input file DEMs.
```matlab
% Fills all depressions on the test DEM, Feldun.tif.
[fill_dem, fill_flow_direction, fill_flow_direction_parents, fill_pits, fill_flow_accumulation] = sdfa(); 

% Fills all depressions on the file someOtherDEM.tif.
[fill_dem, fill_flow_direction, fill_flow_direction_parents, fill_pits, fill_flow_accumulation] = sdfa('./someOtherDEM.tif') 
  
 %Applies 150 mm of rainfall to the file .someOtherDEM.tif.
 [fill_dem, fill_flow_direction, fill_flow_direction_parents, fill_pits, fill_flow_accumulation] = sdfa('./someOtherDEM.tif', 150) 
```
The following variables are provided as output:

 - fill_dem: filled dem
 - fill_flow_direction: flow direction grid corresponding to the filled DEM state
 - fill_flow_direction: flow direction parents cell array corresponding to the filled DEM state
 - fill_pits: pit ID label grid corresponding to the filled DEM state
 - fill_flow_accumulation: flow accumulation grid corresponding to the filled DEM state


#### Important Note:  
**It is recommended that MATLAB be [allocated additional Java Heap Memory](https://www.mathworks.com/help/matlab/matlab_external/java-heap-memory-preferences.html) resources and [parallel workers](https://www.mathworks.com/help/parallel-computing/parallel-preferences.html) to improve computational performance.**


This repo is comprised of the following scripts:

 1. FullAlgorithm.m
 2. d8FlowDirection.m
 3. d8FlowDirectionParents.m
 4. getDepressions.m
 5. fillDepressions.m
 6. flowAccumulation.m


### FullAlgorithm.m
As opposed to running the algorithm as a function with fixed outputs (`sdfa.m`), `FullAlgorithm.m` can be used to get under the hood of the algorithm and explore alternative and other outputs. This is the main script file which calls methods to compute flow direction, initialize the set of depressions, fill depressions sequentially, then generate flow accumulation. This file inside of the `generic/` directory is used to process an unmasked DEM file. By contrast, FullAlgorithm.m scripts within the `manuscript/` directory is used to process a DEM file masked to the shape of the study site watershed (a subset of the Cedar Creek watershed, in this case) for the purposes of recording additional outputs specific to that outlet point which is not located along the border of the DEM file. 

### d8FlowDirection.m
This script computes flow direction using a new implementation whereby flow direction is expressed as the destination index to which that cell flows (using linear cell indexing). This outputs an integer matrix called `flow_direction`. 

### d8FlowDirectionParents.m
By computing the flow direction "parents" of each cell, we generate a data structure which provides for more rapid upstream traversals of flow direction data, essential for computing depression parameters, contributing area, etc. 

### getDepressions.m
This script generates the pit ID label matrix which is the size of the DEM, assigning each cell an integer value equal to the ID of the depression to which it belongs. This also generates the following parameters for each depression:

 - pitId: ID value assigned to this depression in the `pits` and `fill_pits` ID label grids
 - pitCell: pit bottom cell, assigned a flow direction of -1. All cells belonging to this depression flow to this cell.
 - areaCellCount: the integer number of cells comprising this depression
 - spilloverElevation: water level elevation (in meters) in which this depression will overflow into an adjacent depression
 - vca: volume-to-contributing area ratio, equivalent to rainfall excess depth necessary to funnel enough water into this depression to fill it (does not account for infiltration and other losses)
 - volume: the volume of water that this depression retains (cubic meters)
 - filledVolume: zero initially, this value tallies up the volume of water retained over previous iterations of the filling process.
 - cellOverflowInto: cell in an adjoining depression to which this depression will overflow
 - pairs: paired border cells (one just inside the perimeter of this depression, one cell belonging to an adjacent depression) used in `fillDepressions.m` to quickly compute the new spillover location as each depression merges.
 - cellIndexes: array of cell indexes belonging to this depression;
 - 

### fillDepressions.m
This script fills depressions located in the DEM via `getDepressions.m` in an iterative, sequential manner until either A) all depressions are filled or B) the specified rainfall excess amount is applied where remaining depressions are left unfilled. 

### flowAccumulation.m
This script is used to compute flow accumulation.

# Manuscript Results Replication
For manuscript-related details, see [here](https://github.com/WaterApps/SequentialDepressionFillingAlgorithm_Matlab/manuscript).
