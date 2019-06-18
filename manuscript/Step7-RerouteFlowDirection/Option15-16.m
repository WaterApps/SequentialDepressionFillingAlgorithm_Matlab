function[dem, flow_direction, pits, rainfall_excess, runoff_areas, number_of_pits, mean_depth, std_depth, mean_area, std_area, storage_volume, runoff_volume, depthFlow] = fillDepressionsReRouteOptions4Direct(dem, flow_direction, pits, pairs, boundaryCells, cellIndexes, pitId, pitCell, areaCellCount, spilloverElevation, vca, volume, filledVolume, cellOverflowInto, perimeterCell, CedarCreekOutletPoint, R, dem_selection)
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
% shape = false(size(dem));
% nanMins = ones(size(shape)).*realmax('single');
depthFlow = nan(size(dem));
idx = 2;
fprintf('Filling depressions...')
first_pit = order(1);
second_pit = pits(cellOverflowInto(first_pit));
jj = find(~isnan(flow_direction) & flow_direction > 0);
ii = double(flow_direction(~isnan(flow_direction) & flow_direction > 0));
flow_direction_parents = cell(size(flow_direction));
for i = 1 : numel(ii)
    flow_direction_parents{ii(i)} = [flow_direction_parents{ii(i)}, jj(i)];
end
while ~isnan(vca(second_pit)) && ~isinf(vca(first_pit))
    finalOrder(idx) = first_pit;
    rainfall_excess(idx) = vca(first_pit);

%% updateFlowDirectionAndParentsAdjacent
    % reorient filled zone back toward spillover cell, change pit ID to
    % second_pit, and raise DEM
    previousCell = cellOverflowInto(first_pit);
    curCell = perimeterCell(first_pit);
    while flow_direction(curCell) > 0
        nextCell = flow_direction(curCell); % store the old flow direction;
        
        % remove current cell from array of flow direction parent of the
        % next cell downstream
        flow_direction_parents{nextCell}(flow_direction_parents{nextCell} == curCell) = [];

        flow_direction(curCell) = previousCell; % set the new flow direction;
        % add new flow direction parents
        if isempty(find(flow_direction_parents{previousCell} == curCell, 1))
            flow_direction_parents{previousCell} = [flow_direction_parents{previousCell}, curCell];
        end
        
        previousCell = curCell; % set for the next iteration
        curCell = nextCell;
    end
%%    
    % re-ID first pit, raise/fill first pit cells
    lessThanSpillover = cellIndexes{first_pit}(dem(cellIndexes{first_pit}) <= spilloverElevation(first_pit));
    depthFlow(lessThanSpillover) = vca(first_pit); % mark depthFlow with vca
    dem(lessThanSpillover) = spilloverElevation(first_pit); % raise DEM
    
    areaCellCount(second_pit) = areaCellCount(first_pit) + areaCellCount(second_pit);

    % update pits matrix
    pits([cellIndexes{first_pit}]) = second_pit; % update pit matrix
    
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
    % intersection (lines no longer along the perimeter)
    pairs{second_pit} = [pairs{first_pit}; pairs{second_pit}];
    pairs{second_pit} = pairs{second_pit}(any(pits(pairs{second_pit}) ~= second_pit, 2), :);

    % get minimum spillover elevation and cellOverflowInto
    [val, ord] = min(max(dem(pairs{second_pit}(:, 1:2)), [], 2));
    spilloverElevation(second_pit) = val;
    cellOverflowInto(second_pit) = pairs{second_pit}(ord, pits(pairs{second_pit}(ord, :)) ~= second_pit);
    perimeterCell(second_pit) = pairs{second_pit}(ord, pits(pairs{second_pit}(ord, :)) == second_pit);

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
profile viewer
pitId = pitId(finalOrder);
pitCell = pitCell(finalOrder, :);
areaCellCount = areaCellCount(finalOrder);
% spilloverElevation = spilloverElevation(finalOrder);
vca = vca(finalOrder);
volume = volume(finalOrder);
filledVolume = filledVolume(finalOrder);
cellOverflowInto = cellOverflowInto(finalOrder);

function child = recursiveMarkPits(i, oldId, newId, spillElev, depth)
    fd = flow_direction(i);
    if (dem(i) < spillElev)
        dem(i) = spillElev;
        depthFlow(i) = depth;
    end
    pits(i) = newId;
    if (fd == -1) % reached pit bottom
        return
    elseif (pits(i) == oldId) % downstream is yet unidentified
        child = recursiveMarkPits(fd, oldId, newId, spillElev, depth);
    end
    return
end

end