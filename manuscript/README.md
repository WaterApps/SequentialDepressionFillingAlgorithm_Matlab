### Manuscript Title
Optimizations of a sequential depression-filling algorithm and results on a watershed in northeastern Indiana, USA

### Author details:
Name: Samuel A. Noel

Organization: Purdue University Department of Agricultural and Biological Engineering, West Lafayette, IN, USA

Contact: see [my profile](https://github.com/sanoel)


## Manuscript Results Replication
The `/manuscript` directory can be used to replicate results provided in the manuscript. Also included are previously recorded profiles which document the performance results for each respective combination and DEM resolution. To run, open `FullAlgorithm.m` and assign a value of `1`, `2`, or `3` to the `dem_selection` variable corresponding to the 30-meter, 10-meter, and 3-meter DEMs, respectively. Then run the script. Note, depending on the combination and DEM selected, run times may be very long.

As far as test data, the included DEMs `CedarUpper30m.tif`, `CedarUpper10m.tif`, and `CedarUpper3m.tif` are referenced in the `FullAlgorithm.m` scripts as necessary. In the `/generic`, an unmasked DEM `Feldun.tif` is provided as a test DEM.


