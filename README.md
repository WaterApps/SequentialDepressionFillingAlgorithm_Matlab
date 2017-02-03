# Sequential Depression-Filling Algorithm

##  Instructions
1. `FullAlgorithmCedarCreek.m` is the main script. Using a sequential depression-filling algorithm, it generates several outputs for the upper portion of the Cedar Creek Watershed in northeastern Indiana.
2. At the top of the script, under the `%% INPUTS` heading, `dem_selection` should be set to `1`, `2`, or `3` based on whether the user wants to run the algorithm on the 30-meter, 10-meter, or 3-meter-resolution DEM, respectively.
3. Run the script. The 30-meter DEM takes several minutes, the 10-meter several hours, and the 3-meter DEM several days depending on the particular machine.
4. There are several outputs:
a. Matlab profile of the script performance.
b. `.png` images corresponding to state of the `pits` variable at each of the `rainfall_excess_frames` specified in `fillDepressionsCedarCreek.m`. E.g., '30m_15048.png' corresponds to the 30-meter DEM on the 15048th iteration.
c. The full workspace, named `CedarUpper_<30m/10m/03m>_outputs`, containing the DEM before and after filling depressions, flow direction before and after filling depressions, the pit identification matrix before and after filling, as well as several variables storing values at each iteration of the filling process such as rainfall excess, contributing area, etc.
5. To generate several plots using algorithm outputs, run `generatePlots.m`
## License
Apache 2.0

Copyright 2017 Samuel Noel

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
