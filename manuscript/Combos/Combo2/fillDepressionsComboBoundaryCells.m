function[times, dem, flow_direction, pits, rainfall_excess, runoff_areas, number_of_pits, mean_depth, std_depth, mean_area, std_area, storage_volume, runoff_volume, depthFlow] = ...
    fillDepressionsComboBoundaryCells(dem, flow_direction, pits, boundaryCells, pitId, pitCell, areaCellCount, spilloverElevation, vca, volume, filledVolume, cellOverflowInto, CedarCreekOutletPoint, R, dem_selection)
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
[numrows, numcols] = size(dem);
cedarCreekPitId = pits(CedarCreekOutletPoint(1), CedarCreekOutletPoint(2));
vca(cedarCreekPitId) = Inf;
rainfall_excess_frames = [0, 25, 75, 150, 300, 500, NaN]./1000; % meters
rainfall_excess_frames_idx = 2;
cellsize = R.CellExtentInWorldX; % DEM cellsize in meters

% Order pits according to spillover time and initialize the current maximum
% and minimum pit ID numbers
[~, order] = sort(vca);
nonNanCount = nansum(nansum(~isnan(dem)));
disp(strcat('Cells with defined values: ', num2str(nonNanCount)));

% Find the total number of iterations/merges
potential_merges = max(pitId);
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
times = zeros(potential_merges, 1);
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
indexes = reshape(1: numel(flow_direction), size(flow_direction));
% ii = double(flow_direction(~isnan(flow_direction) & flow_direction > 0));
% jj = find(~isnan(flow_direction) & flow_direction > 0);
% flow_direction_parents = cell(size(flow_direction));
% for i = 1 : numel(ii)
%     flow_direction_parents{ii(i)} = [flow_direction_parents{ii(i)}, jj(i)];
% end
while idx <= potential_merges
    iTic = tic();
    finalOrder(idx) = first_pit;
    rainfall_excess(idx) = vca(first_pit);
    
    dem(dem <= spilloverElevation(first_pit) & pits == first_pit) = spilloverElevation(first_pit);
    pits(pits==first_pit) = second_pit; % update pit matrix
    
    % reorient filled zone back toward spillover cell, 
    flow_direction(pitCell(first_pit)) = cellOverflowInto(first_pit);
%     flow_direction_parents{cellOverflowInto(first_pit)} = [flow_direction_parents{cellOverflowInto(first_pit)}, pitCell(first_pit)];
    areaCellCount(second_pit) = areaCellCount(first_pit) + areaCellCount(second_pit);

    if idx == potential_merges
        times(idx-1) = toc(iTic);
        break;
    end
    
    % find the spillover of the new combined depression; prune out those
    % which are no longer on the boundary
    boundaryCells{second_pit} = [boundaryCells{first_pit}; boundaryCells{second_pit}];
    spilloverElevation(second_pit) = realmax('double');
    for p = length(boundaryCells{second_pit}) : -1 : 1 % reverse order so that we can prune in the same pass
        stillBoundary = false;
        [r, c] = ind2sub(size(dem), boundaryCells{second_pit}(p));
        for x = -1 : 1 % loop through neighboring cells
            for y = -1 : 1
                if x == 0 && y ==0
                    continue; % skip center cell of 3x3 neighborhood
                end
                if r+y > numrows || r+y < 1 || c+x > numcols || c+x < 1
                    continue; % skip neighbors outside the dem range
                end

                neighbor = indexes(r+y, c+x);
                curCell = boundaryCells{second_pit}(p);

                % First, check for neighbors to check next.
                if pits(neighbor) ~= pits(curCell)
                    stillBoundary = true;
                    if isnan(pits(neighbor))
                        continue;
                    end
                    elev = max(dem([neighbor, curCell])); % lines are always oriented starting at lowest linear index
                    if (elev < spilloverElevation(second_pit))
                        spilloverElevation(second_pit) = elev;
                        cellOverflowInto(second_pit) = neighbor;
                    end
                end
            end
        end
        if (~stillBoundary)%prune out boundaries no longer on the perimeter
            boundaryCells{second_pit}(p) = [];
        end
    end

    % Now compute volume of the pit
    % compute filled volume and initialize volume calculation (more below)
    filledVolume(second_pit) = volume(first_pit) + filledVolume(second_pit);
    volume(second_pit) = filledVolume(second_pit);
    volume(second_pit) = volume(second_pit) + sum((spilloverElevation(second_pit) - (dem(pits == second_pit & dem <= spilloverElevation(second_pit)))).*cellsize*cellsize);
    
    vca(second_pit) = volume(second_pit)/((cellsize^2).*areaCellCount(second_pit));
    if vca(second_pit) < 0
        vca(second_pit) = Inf;
    end
    
    vca(first_pit) = NaN;
    cedarCreekPitId = pits(CedarCreekOutletPoint(1), CedarCreekOutletPoint(2));
    vca(cedarCreekPitId) = Inf;
    boundaryCells{first_pit} = [];
    
    [~, order] = sort(vca);
    idx = idx + 1
    first_pit = order(1);
    second_pit = pits(cellOverflowInto(first_pit));
    times(idx-1) = toc(iTic);
end
end
