function[dem, flow_direction, pits, rainfall_excess, runoff_areas, number_of_pits, mean_depth, std_depth, mean_area, std_area, storage_volume, runoff_volume, depthFlow] = fillDepressionsVolumeElevArrayNoIndexes(dem, flow_direction, pits, pairs, elevations, cellIndexes, pitId, pitCell, areaCellCount, spilloverElevation, vca, volume, filledVolume, cellOverflowInto, CedarCreekOutletPoint, R, dem_selection)
% Fill depressions in the CedarUpper DEMs. Generate output values and
% images.
%
% Outputs: 
% dem - elevation matrix
%
% flow_direction - flow direction matrix (see d8FlowDirection.m)
%
% flow_direction_parents - flow direction parent matrix (see
% d8FlowDirection.m)
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
rainfall_excess_frames = [0, 25, 75, 150, 300, 500, NaN]./1000; % meters
rainfall_excess_frames_idx = 2;
cellsize = R.CellExtentInWorldX; % DEM cellsize in meters
newtic = tic();
% Order pits according to spillover time and initialize the current maximum
% and minimum pit ID numbers
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
% shape = false(size(dem));
% nanMins = ones(size(shape)).*realmax('single');
depthFlow = nan(size(dem));
idx = 2;
fprintf('Filling depressions...')
first_pit = order(1);
second_pit = pits(cellOverflowInto(first_pit));
 
while ~isnan(vca(second_pit)) && ~isinf(vca(first_pit))
    finalOrder(idx) = first_pit;
    rainfall_excess(idx) = vca(first_pit);

    % re-ID first pit, raise/fill first pit cells
    depthFlow(dem < spilloverElevation(first_pit) & pits == first_pit & isnan(depthFlow)) = vca(first_pit);
    dem(pits == first_pit & dem < spilloverElevation(first_pit)) = spilloverElevation(first_pit);
    areaCellCount(second_pit) = areaCellCount(first_pit) + areaCellCount(second_pit);

    % fix flow direction by pointing pit bottom towards outlet
    flow_direction(pitCell(first_pit)) = cellOverflowInto(first_pit);
    pits(pits == first_pit) = second_pit;
    
    if (idx == length(pitCell))
        number_of_pits(idx) = sum(~isnan(vca));
        runoff_areas(idx) = areaCellCount(cedarCreekPitId).*100./nonNanCount;
        mean_depth(idx) = nanmean(vca);
        std_depth(idx) = nanstd(vca);
        mean_area(idx) = nanmean(areaCellCount)*0.0001*cellsize^2;
        std_area(idx) = nanstd(areaCellCount)*0.0001*cellsize^2;
        storage_volume(idx) = ((nonNanCount - areaCellCount(cedarCreekPitId))*rainfall_excess(idx)*cellsize^2) + filledVolume(cedarCreekPitId);
        runoff_volume(idx) = (areaCellCount(cedarCreekPitId)*rainfall_excess(idx)*cellsize^2) - filledVolume(cedarCreekPitId);
        break 
    end
    
    % Take the union of the sets of boundary lines; prune out the 
    % intersection (those lines are no longer along the perimeter)
    pairs{second_pit} = [pairs{first_pit}; pairs{second_pit}];
    [~, b, c] = unique(pairs{second_pit}(:, 1:2), 'rows'); % get unique values and indexes
    d = accumarray(c, 1); % get counts of each unique value;
    pairs{second_pit} = pairs{second_pit}(b(d<2), :); % for any duplicates (d >= 2), remove both entries 

    [val, ord] = min(max(dem(pairs{second_pit}(:, 1:2)), [], 2)); % get minimum spillover elevation
    spilloverElevation(second_pit) = val;
    cellOverflowInto(second_pit) = pairs{second_pit}(ord, 3);

    filledVolume(second_pit) = volume(first_pit) + filledVolume(second_pit);
    volume(second_pit) = volume(second_pit) + volume(first_pit);
    
    elevations{second_pit} = [elevations{second_pit}; elevations{first_pit}];
    lessThans = elevations{second_pit}(:, 1) <= spilloverElevation(second_pit);
    firstIdx = find(lessThans, 1);
    volume(second_pit) = volume(second_pit) + sum((spilloverElevation(second_pit) - elevations{second_pit}(lessThans, 1)).*elevations{second_pit}(lessThans, 2)*cellsize*cellsize);
    count = sum(elevations{second_pit}(lessThans, 2));
    elevations{second_pit}(lessThans, :) = 0;
    elevations{second_pit}(firstIdx, :) = [spilloverElevation(second_pit), count];    
    elevations{second_pit} = elevations{second_pit}(any(elevations{second_pit}, 2), :);
    lessThans = [];


    vca(second_pit) = volume(second_pit)/((cellsize^2).*areaCellCount(second_pit));
    if vca(second_pit) < 0
        vca(second_pit) = Inf;
    end
  
    vca(first_pit) = NaN;
    cedarCreekPitId = pits(CedarCreekOutletPoint(1), CedarCreekOutletPoint(2));
    vca(cedarCreekPitId) = Inf;
    [~, order] = sort(vca);
    pairs{first_pit} = []; % deallocate memory
    elevations{first_pit} = []; % deallocate memory
    
    number_of_pits(idx) = sum(~isnan(vca));
    runoff_areas(idx) = areaCellCount(cedarCreekPitId).*100./nonNanCount;
    mean_depth(idx) = nanmean(vca);
    std_depth(idx) = nanstd(vca);
    mean_area(idx) = nanmean(areaCellCount)*0.0001*cellsize^2;
    std_area(idx) = nanstd(areaCellCount)*0.0001*cellsize^2;
    storage_volume(idx) = ((nonNanCount - areaCellCount(cedarCreekPitId))*rainfall_excess(idx)*cellsize^2) + filledVolume(cedarCreekPitId);
    runoff_volume(idx) = (areaCellCount(cedarCreekPitId)*rainfall_excess(idx)*cellsize^2) - filledVolume(cedarCreekPitId);

    if (rainfall_excess(idx) >= rainfall_excess_frames(rainfall_excess_frames_idx))
        RGB = ind2rgb(pits, pitsColormap);
        image(RGB);
        axis equal;
        set(gca,'visible','off');
        set(gca,'position',[0 0 1 1], 'units', 'normalized');
        drawnow;
        imwrite(RGB, strcat(dem_selection, num2str(idx),'.png'));
        geotiffwrite(strcat(dem_selection, num2str(idx), '.tif'), pits, R, 'CoordRefSysCode', 26916);
        rainfall_excess_frames_idx = rainfall_excess_frames_idx + 1;
    end
    
    idx = idx + 1;
%     toc(newtic)
    first_pit = order(1);
    second_pit = pits(cellOverflowInto(first_pit));
end

size(pitId)
size(finalOrder)
pitId = pitId(finalOrder);
pitCell = pitCell(finalOrder, :);
areaCellCount = areaCellCount(finalOrder);
% spilloverElevation = spilloverElevation(finalOrder);
vca = vca(finalOrder);
volume = volume(finalOrder);
filledVolume = filledVolume(finalOrder);
cellOverflowInto = cellOverflowInto(finalOrder);
end


