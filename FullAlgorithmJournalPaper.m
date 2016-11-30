clear all;
close all;
profile -memory on;

%% INPUTS
dem_selection = '03m'; % '30m' '10m' '03m'
filename = 'CedarUpper_3m.tif'; %CedarUpper_30m.tif, CedarUpper_10m.tif, CedarUpper_3m.tif
rainfall_depth = 10; % in inches
rainfall_duration = 24; % in hours

%%
if (exist(strcat('CedarUpper_', dem_selection, 'fill_inputs.mat'), 'file') == 2)
    disp('found mat file with filling inputs')
    load(strcat('CedarUpper_', dem_selection, 'fill_inputs.mat'));
else
    disp('filling inputs mat file was not found')
    intensity = (rainfall_depth*0.0254)/rainfall_duration; % meters per hour
    [dem, georef_info] = geotiffread(filename);
    cellsize = georef_info.CellExtentInWorldX; % DEM cellsize in meters
    disp(['DEM loaded...dimensions are: ', num2str(size(dem))]);
    CedarUpperOutletPoints = [685, 470; % outlet point for 30m DEM
                             2052, 1410; % outlet point for 10m DEM
                            6841, 4696]; % outlet point for 3m DEM
    CedarUpperOutletPoint = CedarUpperOutletPoints(dem_selection, :);

    %% Identify NaN values as those with a value of -3.4028e38. Instead, mark as "NaN".
    for i = 1 : numel(dem)
        if dem(i) < 0
            dem(i) = NaN;
        end
    end

    %% Display DEM from data matrix of raw LiDAR Data
    %dem = idwInterpolation(dem); % The DEM has been preprocessed to remove
    % holes and other NaN cells

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
    [pits, pitId, pitCell, areaCellCount, spilloverElevation, spilloverTime, volume, filledVolume, outletCell, cellOverflowInto] = Pits(dem, drainage, flow_direction, flow_direction_parents, cellsize, intensity);
    pittoc = toc(pittic);

    save(strcat('CedarUpper_', dem_selection,'fill_inputs.mat'));
    %% Fill Pits
    disp('Filling Pits')
end
filltic = tic;
[fill_dem, fill_flow_direction, fill_flow_direction_parents, fill_pits, rainfall_excess, runoff_areas, number_of_pits, mean_depth, std_depth, mean_area, std_area, areas, storage_volume, runoff_volume] = fillPitsRecordData(dem, flow_direction, flow_direction_parents, pits, pitId, pitCell, areaCellCount, spilloverElevation, spilloverTime, volume, filledVolume, outletCell, cellOverflowInto, cellsize, intensity, CedarUpperOutletPoint, georef_info, dem_selection);
filltoc = toc(filltic);
save(strcat('CedarUpper_', dem_selection)); %CedarCreek_30m_results 'CedarCreek_30m_sample'

profile off;
profsave(profile('info'),strcat('CedarUpper_', dem_selection, '_profile'));
profile viewer;

figure(2);
percent_running_off = runoff_areas.*100./runoff_areas(end);
stairs(rainfall_excess.*1000, percent_running_off, 'green');
title('Percent of DEM Running Off vs Rainfall Excess');
legend('30m DEM');
ylabel('Percent Area Running Off (%)');
xlabel('Rainfall Excess (mm)');
xlim([0, 500]);
ylim([0, 100]);

% Plot of storage (%)
figure(3);
percent_storage = storage_volume.*100./storage_volume(end);
stairs(rainfall_excess.*1000, percent_storage, '-r', 'LineWidth', 1.5);
title('Depression Storage (Percent of Potential Storage)');
ylabel('Depression Storage (%)');
xlabel('Rainfall Excess (mm)');
legend('30m DEM', '10m DEM', '3m DEM', 'Location', 'southeast');
xlim([0, 500]);
ylim([0, 100]);

% Plot of storage/runoff (%)
% nonNanCount = nansum(nansum(~isnan(dem)));
storage_runoff = storage_volume.*100./(storage_volume + runoff_volume);
figure(4);
% plot(rainfall_excess.*1000, percent_runoff, '-r', 'LineWidth', 1.5);
area(rainfall_excess.*1000, storage_runoff);
title('Stored Depth vs Runoff Depth');
ylabel('Depression Storage and Runoff Depth (mm)');
xlabel('Rainfall Excess (mm)');
% legend('30m DEM', '10m DEM', '3m DEM', 'Location', 'southeast');
xlim([0, 500]);
dim = [.17 .12 .15 .15];
str = '% Depression Storage';
annotation('textbox',dim,'String',str,'FitBoxToText','on', 'EdgeColor', 'none', 'BackgroundColor', 'none', 'Color', 'white', 'FontWeight', 'bold');
dim = [.50 .5 .15 .15];
str = '% Runoff';
annotation('textbox',dim,'String',str,'FitBoxToText','on', 'EdgeColor', 'none', 'BackgroundColor', 'none', 'FontWeight', 'bold');
% ylim([0, 100]);

