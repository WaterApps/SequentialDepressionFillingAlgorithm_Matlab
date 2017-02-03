clear all;
close all;
profile -memory on;

dem_options = ['30m'; '10m'; '03m']; %correspond to files CedarUpper_30m.tif, CedarUpper_10m.tif, CedarUpper_03m.tif
CedarUpperOutletPoints = [685, 470; % outlet point for 30m DEM
                          2052, 1410; % outlet point for 10m DEM
                          6841, 4696]; % outlet point for 3m DEM
%% INPUTS
dem_selection = 1; % choose DEM from above list "dem_options"

%%
CedarUpperOutletPoint = CedarUpperOutletPoints(dem_selection, :);
if (exist(strcat('CedarUpper_', dem_options(dem_selection, :),'fill_inputs.mat'), 'file') == 2)
    disp(strcat('Filling inputs file CedarUpper_', dem_options(dem_selection, :),'fill_inputs.mat was found'))
    load(strcat('CedarUpper_', dem_options(dem_selection, :),'fill_inputs.mat'));
else
    disp(strcat('Filling inputs file CedarUpper_', dem_options(dem_selection, :),'fill_inputs.mat was not found'))
    [dem, georef_info] = geotiffread(strcat('CedarUpper_', dem_options(dem_selection, :)));
    cellsize = georef_info.CellExtentInWorldX; % DEM cellsize in meters
    disp(['DEM loaded...dimensions are: ', num2str(size(dem))]);

    %% The clipped DEMs have an assigned NaN value of -3.4028e38. Find and mark as NaN.
    for i = 1 : numel(dem)
        if dem(i) < -100
            dem(i) = NaN;
        end
    end

    %% Compute Flow Direction Matrix/Map
    disp('Computing Flow Direction');
    flowtic = tic;
    [flow_direction, flow_direction_parents] = d8FlowDirection(dem);
    flowtoc = toc(flowtic);

    %% Identify Pits, Compute Matrix/Map with Pit ID for each cell
    disp('Creating Pit Database');
    pittic = tic;
    [pits, pitId, pitCell, areaCellCount, spilloverElevation, vca, volume, filledVolume, outletCell, cellOverflowInto] = getDepressions(dem, flow_direction, flow_direction_parents, georef_info.CellExtentInWorldX);
    pittoc = toc(pittic);

    save(strcat('CedarUpper_', dem_options(dem_selection, :),'fill_inputs.mat'));
end
%% Fill Pits
disp('Filling Pits')
filltic = tic;
[fill_dem, fill_flow_direction, fill_flow_direction_parents, fill_pits, rainfall_excess, runoff_areas, number_of_pits, mean_depth, std_depth, mean_area, std_area, storage_volume, runoff_volume] = fillDepressionsCedarCreek(dem, flow_direction, flow_direction_parents, pits, pitId, pitCell, areaCellCount, spilloverElevation, vca, volume, filledVolume, outletCell, cellOverflowInto, CedarUpperOutletPoint, georef_info, dem_options(dem_selection, :));
filltoc = toc(filltic);

% Save current workspace to save algorithm results.
save(strcat('CedarUpper_', dem_options(dem_selection, :),'_outputs'));

% Save the profile results to analyze algorithm performance
profile off;
profsave(profile('info'),strcat('CedarUpper_', dem_options(dem_selection, :), '_profile'));
profile viewer;


