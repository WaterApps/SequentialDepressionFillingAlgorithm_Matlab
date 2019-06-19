clear all;
close all;
profile -memory on;
try 
%% INPUTS
dem_selection = 1; % choose DEM from above list "dem_options"
dem_options = ['30m'; '10m'; '03m']; %correspond to files CedarUpper_30m.tif, CedarUpper_10m.tif, CedarUpper_03m.tif
CedarUpperOutletPoints = [685, 470; % outlet point for 30m DEM
                          2052, 1410; % outlet point for 10m DEM
                          6841, 4696]; % outlet point for 3m DEM

if (ispc)

  addpath('..\..\..\')
else
  addpath('../../../')
end

%%
CedarUpperOutletPoint = CedarUpperOutletPoints(dem_selection, :);
[dem, georef_info] = geotiffread(strcat('CedarUpper_', dem_options(dem_selection, :)));
cellsize = georef_info.CellExtentInWorldX; % DEM cellsize in meters
disp(['DEM loaded...dimensions are: ', num2str(size(dem))]);

%% The clipped DEMs have an assigned NaN value of -3.4028e38. Find and mark as NaN.
%dem(dem < -100) = NaN;
dem = double(dem);
parfor i = 1 : numel(dem)
    if (dem(i) < -100)
        dem(i) = NaN
    end
end

%% Recompute Flow Direction
flowDirectionFile = strcat('ComboCacheCedarUpper_', dem_options(dem_selection, :), '_FlowDirection.mat');
if (false) %exist(flowDirectionFile, 'file') == 2)
    disp('Flow direction matrix was found')
    load(flowDirectionFile);
else
    disp('Flow direction matrix was NOT found')
    disp('Computing Flow Direction');
    flowtic = tic;
    [flow_direction] = d8FlowDirection(dem);
    flowtoc = toc(flowtic);
    save(flowDirectionFile, 'flow_direction', 'flowtoc')
end
 [flow_direction_parents] = d8FlowDirectionParents(flow_direction);

%% Identify Pits, Compute Matrix/Map with Pit ID for each cell
pitsFile = strcat('ComboCacheCedarUpper_', dem_options(dem_selection, :), '_Pits.mat');
if (false) %exist(pitsFile, 'file') == 2)
    disp('Pits data was found')
    load(pitsFile);
else
    disp('Creating Pit Database');
    [pits, pitId, pitCell, areaCellCount, spilloverElevation, vca, volume, filledVolume, cellOverflowInto, pairs, cellIndexes] = getDepressions2(dem, flow_direction, flow_direction_parents, georef_info.CellExtentInWorldX);
    save(pitsFile, 'pits', 'pitId', 'pitCell', 'areaCellCount', 'spilloverElevation', 'vca', 'volume', 'filledVolume', 'cellOverflowInto', 'pairs', 'cellIndexes')
end

%% Fill Pits
disp('Filling Pits')
fillFile = strcat('ComboCacheCedarUpper_', dem_options(dem_selection, :), '_outputs.mat');
if (false) %exist(fillFile, 'file') == 2)
    disp('Fill data was found')
    load(fillFile);
else
    disp('Creating Fill Outputs');
    filltic = tic;
    [times, fill_dem, fill_flow_direction, fill_pits, rainfall_excess, runoff_areas, number_of_pits, mean_depth, std_depth, mean_area, std_area, storage_volume, runoff_volume, depthFlow, first_pit_areas, second_pit_areas] = ...
        fillDepressionsComboCache(dem, flow_direction, pits, pairs, cellIndexes, pitId, pitCell, areaCellCount, spilloverElevation, vca, volume, filledVolume, cellOverflowInto, CedarUpperOutletPoint, georef_info, dem_options(dem_selection, :));
    filltoc = toc(filltic);
    save(fillFile);
end

%% Flow Accumulation
disp('Computing Flow Accumulation')
[fill_flow_accumulation] = flowAccumulation(fill_flow_direction);

% Save the profile results to analyze algorithm performance
profile off;
%profsave(profile('info'),strcat('ComboCacheCedarUpper_', dem_options(dem_selection, :), '_profile'));
profile viewer;
catch ME
    profile off
    rethrow(ME)
end

