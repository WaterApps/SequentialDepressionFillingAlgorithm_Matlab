function[dem, flow_direction, pits, rainfall_excess, runoff_areas, number_of_pits, mean_depth, std_depth, mean_area, std_area, storage_volume, runoff_volume, depthFlow, first_pit_areas, second_pit_areas] = ...
    fillDepressionsSortOrder(dem, flow_direction, pits, pairs, cellIndexes, pitId, pitCell, areaCellCount, spilloverElevation, vca, volume, filledVolume, cellOverflowInto, CedarCreekOutletPoint, R, dem_selection)
% Fill depressions in the CedarUpper DEMs. Generate output values and
% images.
%
% Outputs: 
% dem - elevation matrix
%
% flow_direction - flow direction matrix (see d8FlowDirection.m)
%
% pits - matrix identifying depressions (see getDepressions.m)
%
% rainfall_excess - The rainfall_excess value of each iteration (meters).
%
% runoff_areas - The contributing area of the Cedar Creek watershed at each
% iteration (hectares).
%
% number_of_pits - Number of depressions remaining at each iteration.
%
% mean_depth - Mean vca across all remaining depressions (meters).
%
% std_depth - Standard deviation of the vca values of remaining
% depressions (meters).
%
% mean_area - Mean area (hectares) of remaining depressions at each
% iteration.
%
% std_area - Standard deviation of the areas of remaining depressions
% (hectares)
%
% storage_volume - The total volume stored at each iteration (m^3).
%
% runoff_volume - The total volume that has run off at each iteration (m^3).

%First, find the pit containing the Cedar Creek outlet, and set its
%spillover time to infinity.  This will prevent it from overflowing into
%other depressions. Only other depressions may overflow into it.
cedarCreekPitId = pits(CedarCreekOutletPoint(1), CedarCreekOutletPoint(2));
vca(cedarCreekPitId) = Inf;
% eventualSize = nnz(~isnan(dem));
% currentSize = length(cellIndexes{cedarCreekPitId});
% cellIndexes{cedarCreekPitId}(currentSize+1:eventualSize) = 0; % this thing will eventually hold all of the cell indexes

rainfall_excess_frames = [0, 25, 75, 150, 300, 500, NaN]./1000; % meters
rainfall_excess_frames_idx = 2;
cellsize = R.CellExtentInWorldX; % DEM cellsize in meters

[~, order] = sort(vca);
% pitId = pitId(order);
% pitCell = pitCell(order, :);
% areaCellCount = areaCellCount(order);
% spilloverElevation = spilloverElevation(order);
% vca = vca(order);
% volume = volume(order);
% filledVolume = filledVolume(order);
% cellOverflowInto = cellOverflowInto(order);
max_id = max(pitId);
nonNanCount = nansum(nansum(~isnan(dem)));
disp(strcat('Cells with defined values: ', num2str(nonNanCount)));

% Find the total number of iterations/merges
potential_merges = max_id;
% Preallocate arrays of variables that will be used in analysis
first_pit_areas = zeros(potential_merges, 1);
second_pit_areas = zeros(potential_merges, 1);
rainfall_excess = zeros(potential_merges, 1); % meters
runoff_areas = zeros(potential_merges, 1); % hectares
number_of_pits = zeros(potential_merges, 1); 
mean_depth = zeros(potential_merges, 1); % meters
std_depth = zeros(potential_merges, 1); % meters
mean_area = zeros(potential_merges, 1); % hectares
std_area = zeros(potential_merges, 1); % hectares
storage_volume = zeros(potential_merges, 1); % cubic meters
runoff_volume = zeros(potential_merges, 1); % cubic meters
finalOrder = nan(potential_merges, 1);

% Get initial values
rainfall_excess(1) = 0;
number_of_pits(1) = sum(~isnan(vca));
runoff_areas(1) = areaCellCount(cedarCreekPitId).*100./nonNanCount;
mean_depth(1) = nanmean(vca);
std_depth(1) = nanstd(vca);
mean_area(1) = nanmean(areaCellCount)*0.0001*cellsize^2;
std_area(1) = nanstd(areaCellCount)*0.0001*cellsize^2;
runoff_volume(1) = 0;
storage_volume(1) = 0;
finalOrder(1) = order(1);

pitsColormap = rand(max(max(pits))+1, 3); % random color for every depression, plus one for depression ID 0.
pitsColormap(1,:) = 1; % make pit ID 0 white
RGB = ind2rgb(pits, pitsColormap);
image(RGB);
axis equal;
set(gca,'visible','off');
set(gca,'position',[0 0 1 1], 'units', 'normalized');
drawnow;
imwrite(RGB, strcat(dem_selection,num2str(1),'.png'));
geotiffwrite(strcat(dem_selection, num2str(1), '.tif'), pits, R, 'CoordRefSysCode', 26916);

depthFlow = nan(size(dem));
idx = 2;
fprintf('Filling depressions...')
first_pit = order(1);
second_pit = pits(cellOverflowInto(first_pit));
while idx <= potential_merges
    finalOrder(idx) = first_pit;
    rainfall_excess(idx) = vca(first_pit);
    first_pit_areas(idx-1) = areaCellCount(first_pit);
    second_pit_areas(idx-1) = areaCellCount(second_pit);
    
    % re-ID first pit, raise/fill first pit cells
    lessThanSpillover = cellIndexes{first_pit}(dem(cellIndexes{first_pit}) <= spilloverElevation(first_pit));
    depthFlow(lessThanSpillover) = vca(first_pit); % mark depthFlow with vca
    dem(lessThanSpillover) = spilloverElevation(first_pit); % raise DEM

    %update flow directions
    flow_direction(pitCell(first_pit)) = cellOverflowInto(first_pit);

    areaCellCount(second_pit) = areaCellCount(first_pit) + areaCellCount(second_pit);
    pits([cellIndexes{first_pit}]) = second_pit; % update pit matrix
   
    % Find spillover elevation and location; append pairs, prune, and find
    pairs{second_pit} = [pairs{first_pit}; pairs{second_pit}];
    pairs{second_pit} = pairs{second_pit}(any(pits(pairs{second_pit}) ~= second_pit, 2), :);
    [val, ord] = min(max(dem(pairs{second_pit}(:, 1:2)), [], 2));
    
    if idx == potential_merges
        break;
    end
    
    spilloverElevation(second_pit) = val;
    cellOverflowInto(second_pit) = pairs{second_pit}(ord, pits(pairs{second_pit}(ord, :)) ~= second_pit);

    % compute filled volume and initialize volume calculation (more below)
    filledVolume(second_pit) = volume(first_pit) + filledVolume(second_pit);
    volume(second_pit) = filledVolume(second_pit);
    
    cellIndexes{second_pit} = [cellIndexes{second_pit}, cellIndexes{first_pit}];
    lessThans = cellIndexes{second_pit}(dem(cellIndexes{second_pit}) <= spilloverElevation(second_pit));

    volume(second_pit) = volume(second_pit) + sum((spilloverElevation(second_pit) - dem(lessThans) )*cellsize*cellsize);
    lessThans = [];

    vca(second_pit) = volume(second_pit)/((cellsize^2).*areaCellCount(second_pit));
    if vca(second_pit) < 0
        vca(second_pit) = Inf;
    end
  
    vca(first_pit) = NaN;
    cedarCreekPitId = pits(CedarCreekOutletPoint(1), CedarCreekOutletPoint(2));
    vca(cedarCreekPitId) = Inf;
    
    [~, order] = sort(vca);
    cellIndexes{first_pit} = [];
    pairs{first_pit} = [];
    
    idx = idx + 1;
    first_pit = order(1);
    second_pit = pits(cellOverflowInto(first_pit));
end

profile viewer
pitId = pitId(finalOrder);
pitCell = pitCell(finalOrder, :);
areaCellCount = areaCellCount(finalOrder);
%spilloverElevation = spilloverElevation(finalOrder);
vca = vca(finalOrder);
volume = volume(finalOrder);
filledVolume = filledVolume(finalOrder);
cellOverflowInto = cellOverflowInto(finalOrder);
end


