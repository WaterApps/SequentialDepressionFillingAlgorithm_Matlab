function[dem, flow_direction, pits] = ...
    fillDepressionsVolume(dem, flow_direction, pits, pairs, elevations, cellIndexes, pitId, pitCell, areaCellCount, spilloverElevation, vca, volume, filledVolume, cellOverflowInto, perimeterCell, CedarCreekOutletPoint, R, dem_selection)
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
cellsize = R.CellExtentInWorldX; % DEM cellsize in meters

max_id = max(pitId);
nonNanCount = nansum(nansum(~isnan(dem)));
disp(strcat('Cells with defined values: ', num2str(nonNanCount)));

% Find the total number of iterations/merges
potential_merges = max_id;

idx = 2;
fprintf('Filling depressions...')
[~, first_pit] = min(vca);
second_pit = pits(cellOverflowInto(first_pit));
ii = double(flow_direction(~isnan(flow_direction) & flow_direction > 0));
jj = find(~isnan(flow_direction) & flow_direction > 0);
flow_direction_parents = cell(size(flow_direction));
for i = 1 : numel(ii)
    flow_direction_parents{ii(i)} = [flow_direction_parents{ii(i)}, jj(i)];
end
[numrows, numcols] = size(dem);
indexes = reshape(1: numel(flow_direction), size(flow_direction));
while idx <= potential_merges    
    previousCell = cellOverflowInto(first_pit);
    curCell = perimeterCell(first_pit);
    while ~isnan(curCell)
        nextCell = flow_direction(curCell); % store the old flow direction;
        flow_direction(curCell) = previousCell; % set the new flow direction;
        
        % add new flow direction parents
        if isempty(find(flow_direction_parents{previousCell} == curCell, 1))
            flow_direction_parents{previousCell} = [flow_direction_parents{previousCell}, curCell];
        end
        
        % remove current cell from array of flow direction parent of the
        % next cell downstream
        if nextCell ~= -1
            flow_direction_parents{nextCell}(flow_direction_parents{nextCell} == curCell) = [];
        else
            break;
        end
           
        previousCell = curCell; % set for the next iteration
        curCell = nextCell;
    end
    
    % re-ID first pit, raise/fill first pit cells
    lessThanSpillover = cellIndexes{first_pit}(dem(cellIndexes{first_pit}) <= spilloverElevation(first_pit));
    dem(lessThanSpillover) = spilloverElevation(first_pit); % raise DEM
    areaCellCount(second_pit) = areaCellCount(first_pit) + areaCellCount(second_pit);

    % update flow direction by pointing pit bottom towards outlet
%     flow_direction(pitCell(first_pit)) = cellOverflowInto(first_pit);
%     flow_direction_parents{cellOverflowInto(first_pit)} = [flow_direction_parents{cellOverflowInto(first_pit)}, pitCell(first_pit)];
    pits([cellIndexes{first_pit}]) = second_pit; % update pit matrix
    
    if idx == potential_merges
        break;
    end
    
    % Take the union of the sets of boundary lines; prune out the 
    % intersection (lines no longer along the perimeter)
    cellIndexes{second_pit} = [cellIndexes{second_pit}, cellIndexes{first_pit}];
    pairs{second_pit} = [pairs{first_pit}; pairs{second_pit}];
    pairs{second_pit} = pairs{second_pit}(any(pits(pairs{second_pit}) ~= second_pit, 2), :);
    [val, ord] = min(max(dem(pairs{second_pit}(:, 1:2)), [], 2));


    spilloverElevation(second_pit) = val;
    cellOverflowInto(second_pit) = pairs{second_pit}(ord, pits(pairs{second_pit}(ord, :)) ~= second_pit);
    perimeterCell(second_pit) = pairs{second_pit}(ord, pits(pairs{second_pit}(ord, :)) == second_pit);

    % compute filled volume and initialize volume calculation (more below)
    filledVolume(second_pit) = volume(first_pit) + filledVolume(second_pit);
    volume(second_pit) = filledVolume(second_pit);
    indicesToCheck = zeros(areaCellCount(second_pit), 1);
    indicesToCheck(1) = pitCell(second_pit);
    j = 1;
    i = 1;
    volume(second_pit) = 0;
    while i <= j
        volume(second_pit) = volume(second_pit) + ((spilloverElevation(second_pit) - double(dem(indicesToCheck(i))))*cellsize*cellsize);
        if dem(indicesToCheck(i)) >= spilloverElevation(second_pit)
            break
        end
        parents = flow_direction_parents{indicesToCheck(i)};
        lessThans = parents(dem(parents) < spilloverElevation(second_pit));
        k = sum(length(lessThans));
        indicesToCheck(j+1:j+k) = lessThans;
        j = j + k;
        i = i + 1;
    end

    
    vca(second_pit) = volume(second_pit)/((cellsize^2).*areaCellCount(second_pit));
    if vca(second_pit) < 0
        vca(second_pit) = Inf;
    end
  
    vca(first_pit) = NaN;
    cedarCreekPitId = pits(CedarCreekOutletPoint(1), CedarCreekOutletPoint(2));
    vca(cedarCreekPitId) = Inf;
    cellIndexes{first_pit} = [];
    
    idx = idx + 1;
    [~, first_pit] = min(vca);
    second_pit = pits(cellOverflowInto(first_pit));
end

[fill_flow_accumulation] = flowAccumulation(flow_direction);
fprintf('Done')
fprintf('\n')

end


