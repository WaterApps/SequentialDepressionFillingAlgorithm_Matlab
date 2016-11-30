clear all;
close all;
%profile -memory on;

%% INPUTS
rainfall_depth = 10; % in inches
rainfall_duration = 24; % in hours
filename = 'CedarUpper_10m.tif'; %CedarUpper_30m.tif, CedarUpper_10m.tif, CedarUpper_3m.tif
dem_selection = 2;
%%
intensity = (rainfall_depth*0.0254)/rainfall_duration; % convert depth to meters, and calculate intensity in meters per hour
[dem, georef_info] = geotiffread(filename);
cellsize = georef_info.CellExtentInWorldX; % DEM cellsize in meters
disp(['DEM loaded...dimensions are: ', num2str(size(dem))]);

% dem = dem(1:200, size(dem, 2)-200:size(dem,2));
% a = nan(300,300);
% a(50:249, 50:250) = dem;
% dem = a;
% imagery =
% imagery = imagery(1:30*200, size(dem, 2)*30-(200*30:size(dem,2)*30));
% 249, 109; % outlet point for sample DEM
CedarUpperOutletPoints = [685, 470; % outlet point for 30m DEM
                         2052, 1410; % outlet point for 10m DEM
                        6733, 4056]; % outlet point for 3m DEM
CedarUpperOutletPoint = CedarUpperOutletPoints(dem_selection, :);

% CedarCreekOutletPoint = [249,109]; % outlet point for sample DEM
% CedarCreekOutletPoint = [1236, 862]; % outlet point for 30m DEM
% CedarCreekOutletPoint = [3695, 2584]; % outlet point for 10m DEM
% CedarCreekOutletPoint = [1236, 862]; % outlet point for 1.5m DEM

%% Identify NaN values as those with a value of -3.4028e38. Instead, mark as "NaN".
for i = 1 : numel(dem)
    if dem(i) < -100
        dem(i) = NaN;
    end
end

%% Display DEM from data matrix of raw LiDAR Data
%dem = idwInterpolation(dem);

%% Get Drainage Data
drainage = sparse(zeros(size(dem)));
% drainage = placeDrainageFeatures(drainage);

%% Compute Flow Direction Matrix/Map
disp('Computing Flow Direction');
flowtic = tic;
[flow_direction, flow_direction_parents] = d8FlowDirectionDrainage(dem, drainage, intensity);
flowtoc = toc(flowtic);

%% Identify Pits, Compute Matrix/Map with Pit ID for each cell
disp('Creating Pit Database');
pittic = tic;
[pits, pitId, pitCell, areaCellCount, spilloverElevation, spilloverTime, volume, filledVolume, outletCell, cellOverflowInto] = PitsGetImages(dem, drainage, flow_direction, flow_direction_parents, cellsize, intensity);
pittoc = toc(pittic);

%% Fill Pits
disp('Filling Pits')
[fill_dem, fill_flow_direction, fill_flow_direction_parents, fill_pits, fill_rainfall_amounts, runoff_areas, number_of_pits, mean_depth, std_depth, mean_area] = fillPitsRecordDataGetImages(dem, flow_direction, flow_direction_parents, pits, pitId, pitCell, areaCellCount, spilloverElevation, spilloverTime, volume, filledVolume, outletCell, cellOverflowInto, cellsize, intensity, CedarUpperOutletPoint, georef_info, dem_selection);