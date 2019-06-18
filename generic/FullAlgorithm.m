clear all;
close all;
profile -memory on;
try 
    dem_filepath = './Feldun.tif';
    visualize_merging = true; % slows the filling process

    %%
    [dem, georef_info] = geotiffread(dem_filepath);
    cellsize = georef_info.CellExtentInWorldX; % DEM cellsize in meters
    disp(['DEM loaded...dimensions are: ', num2str(size(dem))]);

    %% The clipped DEMs have an assigned NaN value of -3.4028e38. Find and mark as matlab NaN.
    dem = double(dem);
    parfor i = 1 : numel(dem)
        if (dem(i) < -100)
            dem(i) = NaN
        end
    end

    %% Recompute Flow Direction
    disp('Flow direction matrix was NOT found')
    disp('Computing Flow Direction');
    [flow_direction] = d8FlowDirection(dem);
    [flow_direction_parents] = d8FlowDirectionParents(flow_direction);

    %% Identify Pits, Compute Matrix/Map with Pit ID for each cell
    disp('Creating Pit Dataset');
    [pits, pairs, cellIndexes, pitId, pitCell, areaCellCount, spilloverElevation, vca, volume, filledVolume, cellOverflowInto] = ...
        getDepressions(dem, flow_direction, flow_direction_parents, georef_info.CellExtentInWorldX);

    %% Fill Pits
    [times, fill_dem, fill_flow_direction, fill_pits, rainfall_excess, number_of_pits, depthFlow] = ...
        fillDepressions(dem, flow_direction, pits, pairs, cellIndexes, pitId, pitCell, areaCellCount, spilloverElevation, vca, volume, filledVolume, cellOverflowInto, georef_info, visualize_merging);

    %% Flow Accumulation
    disp('Computing Flow Accumulation')
    [fill_flow_accumulation] = flowAccumulation(fill_flow_direction);

    % Save the profile results to analyze algorithm performance
    profile off;
    profile viewer;
catch ME
    profile off
    rethrow(ME)
end

