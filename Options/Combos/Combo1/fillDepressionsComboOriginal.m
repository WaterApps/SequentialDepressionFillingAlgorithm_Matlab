function[times, dem, flow_direction, pits, rainfall_excess, runoff_areas, number_of_pits, mean_depth, std_depth, mean_area, std_area, storage_volume, runoff_volume, depthFlow, first_pit_areas, second_pit_areas] = ...
    fillDepressionsComboOriginal(dem, flow_direction, pits, pitId, pitCell, areaCellCount, spilloverElevation, vca, volume, filledVolume, cellOverflowInto, CedarCreekOutletPoint, R, dem_selection)
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

[~, order] = sort(vca);

rainfall_excess_frames = [0, 25, 75, 150, 300, 500, NaN]./1000; % meters
rainfall_excess_frames_idx = 2;
cellsize = R.CellExtentInWorldX; % DEM cellsize in meters

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
[numrows, numcols] = size(dem);
indexes = reshape(1: numel(flow_direction), size(flow_direction));
while idx <= potential_merges
    iTic = tic();
    finalOrder(idx-1) = first_pit;
    rainfall_excess(idx) = vca(first_pit);
    first_pit_areas(idx-1) = areaCellCount(first_pit);
    second_pit_areas(idx-1) = areaCellCount(second_pit);
    
    %% updateFlowDirectionAndParentsAdjacent
    % reorient filled zone back toward spillover cell, change pit ID to
    % second_pit, and raise DEM
    visited = false(numrows, numcols);
    indicesToCheck = zeros(areaCellCount(first_pit)+1, 2);
    [r, c] = ind2sub(size(flow_direction), cellOverflowInto(first_pit));
    indicesToCheck(1, :) = [r, c];
    j = 1; % overwrite the first element; we start at a cell in the neighboring depression to which it'll overflow
    i = 1;
    while i <= j
        r = indicesToCheck(i, 1);
        c = indicesToCheck(i, 2);
        i = i + 1;
        for x = -1 : 1 % loop through neighboring cells
            for y = -1 : 1
                if x == 0 && y ==0
                    continue; % skip center cell of 3x3 neighborhood
                end
                if r+y > numrows || r+y < 1 || c+x > numcols || c+x < 1
                    continue; % skip neighbors outside the dem range
                end
                neighbor = indexes(r+y, c+x);
                curCell = indexes(r, c);
                if pits(neighbor) == (first_pit) && dem(neighbor) <= spilloverElevation(first_pit) && ~visited(neighbor)
                    visited(neighbor) = true;
                    j = j + 1;
                    indicesToCheck(j, :) = [r+y, c+x];
                    % remove old flow direction parents
                    flow_direction(neighbor) = curCell;
                end
            end
        end
    end
    
    % re-ID first pit, raise/fill first pit cells
    indicesToCheck = zeros(areaCellCount(first_pit), 2);
    [r, c] = ind2sub(size(dem), pitCell(first_pit));
    indicesToCheck(1, :) = [r, c];
    dem(r,c) = spilloverElevation(first_pit);
    pits(r,c) = second_pit;
    j = 1; % overwrite the first element; we start at a cell in the neighboring depression to which it'll overflow
    i = 1;
    [numrows, numcols] = size(dem);
    while i <= j
        r = indicesToCheck(i, 1);
        c = indicesToCheck(i, 2);
        i = i + 1;
        for x = -1 : 1 % loop through neighboring cells
            for y = -1 : 1
                if x == 0 && y ==0
                    continue; % skip center cell of 3x3 neighborhood
                end
                if r+y > numrows || r+y < 1 || c+x > numcols || c+x < 1
                    continue; % skip neighbors outside the dem range
                end
                neighbor = indexes(r+y, c+x);
                if pits(neighbor) == first_pit
                    pits(neighbor) = second_pit;
                    j = j + 1;
                    indicesToCheck(j, :) = [r+y, c+x];
                    
                    if dem(neighbor) <= spilloverElevation(first_pit)
                        dem(neighbor) = spilloverElevation(first_pit);
                    end
                end
            end
        end
    end

    areaCellCount(second_pit) = areaCellCount(first_pit) + areaCellCount(second_pit);
        
    if idx == potential_merges
        times(idx-1) = toc(iTic);
        break;
    end
    %% Find spillover elevation and location
    spilloverElevation(second_pit) = realmax('double');
    indicesToCheck = zeros(areaCellCount(second_pit), 2);
    [r, c] = ind2sub(size(dem), pitCell(second_pit));
    indicesToCheck(1, :) = [r, c];
    j = 1;
    i = 1;
    [numrows, numcols] = size(dem);
    while i <= j
        r = indicesToCheck(i, 1);
        c = indicesToCheck(i, 2);
        i = i + 1;
        for x = -1 : 1 % loop through neighboring cells
            for y = -1 : 1
                if x == 0 && y ==0
                    continue; % skip center cell of 3x3 neighborhood
                end
                if r+y > numrows || r+y < 1 || c+x > numcols || c+x < 1
                    continue; % skip neighbors outside the dem range
                end
                neighbor = indexes(r+y, c+x);
                curCell = indexes(r, c);

                if flow_direction(neighbor) == curCell % by definition has same pit ID
                	j = j + 1;
                    indicesToCheck(j, :) = [r+y, c+x];
                elseif pits(neighbor) ~= pits(curCell) && ~isnan(dem(neighbor))
                    elev = max(dem([neighbor, curCell]));
                    if elev < spilloverElevation(second_pit)
                        spilloverElevation(second_pit) = elev;
                        cellOverflowInto(second_pit) = neighbor;
                    end
                end
            end
        end
    end
    
    %% compute filled volume and initialize volume calculation (more below)
    filledVolume(second_pit) = volume(first_pit) + filledVolume(second_pit);
    volume(second_pit) = filledVolume(second_pit);
    indicesToCheck = zeros(areaCellCount(second_pit), 2);
    [r, c] = ind2sub(size(flow_direction), pitCell(second_pit));
    indicesToCheck(1, :) = [r, c];
    j = 2;
    i = 1;
    while i < j
        r = indicesToCheck(i, 1);
        c = indicesToCheck(i, 2);
        i = i + 1;
        volume(second_pit) = volume(second_pit) + ((spilloverElevation(second_pit) - dem(r,c))*cellsize*cellsize);
        for x = -1 : 1 % loop through neighboring cells
            for y = -1 : 1
                if x == 0 && y ==0
                    continue; % skip center cell of 3x3 neighborhood
                end
                if r+y > numrows || r+y < 1 || c+x > numcols || c+x < 1
                    continue; % skip neighbors outside the dem range
                end

                curCell = indexes(r, c);
                if flow_direction(r+y, c+x) == curCell && dem(r+y, c+x) <= spilloverElevation(second_pit)
                    indicesToCheck(j, :) = [r+y, c+x];
                    j = j + 1;
                end
            end
        end
    end

    vca(second_pit) = volume(second_pit)/((cellsize^2).*areaCellCount(second_pit));
    if vca(second_pit) < 0
        vca(second_pit) = Inf;
    end
  
    vca(first_pit) = NaN;
    cedarCreekPitId = pits(CedarCreekOutletPoint(1), CedarCreekOutletPoint(2));
    vca(cedarCreekPitId) = Inf;
         
    [~, order] = sort(vca);
    idx = idx + 1
    first_pit = order(1);
    second_pit = pits(cellOverflowInto(first_pit));
    times(idx-1) = toc(iTic);
end
profile viewer

fprintf('Done')
fprintf('\n')

end


