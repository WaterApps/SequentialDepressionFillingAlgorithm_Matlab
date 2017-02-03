function[dem, flow_direction, flow_direction_parents, pits, rainfall_excess, runoff_areas, number_of_pits, mean_depth, std_depth, mean_area, std_area, storage_volume, runoff_volume] = fillDepressionsCedarCreek(dem, flow_direction, flow_direction_parents, pits, pitId, pitCell, areaCellCount, spilloverElevation, vca, volume, filledVolume, outletCell, cellOverflowInto, CedarCreekOutletPoint, R, dem_selection)
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
[cedarCreekPitIndex, ~] = find(pitId == cedarCreekPitId);
vca(cedarCreekPitIndex) = Inf;
rainfall_excess_frames = [0, 25, 75, 150, 300, 500, NaN]./1000; % meters
rainfall_excess_frames_idx = 2;
cellsize = R.CellExtentInWorldX; % DEM cellsize in meters

% Order pits according to spillover time and initialize the current maximum
% and minimum pit ID numbers
[~, order] = sort(vca);
pitId = pitId(order);
pitCell = pitCell(order, :);
areaCellCount = areaCellCount(order);
spilloverElevation = spilloverElevation(order);
vca = vca(order);
volume = volume(order);
filledVolume = filledVolume(order);
outletCell = outletCell(order, :);
cellOverflowInto = cellOverflowInto(order, :);
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

% Get initial values
rainfall_excess(1) = 0;
number_of_pits(1) = sum(~isnan(vca));
runoff_areas(1) = areaCellCount(cedarCreekPitIndex).*100./nonNanCount;
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

idx = 2;
while ~isnan(vca(2))
    rainfall_excess(idx) = vca(1);
    second_pit_ID = pits(cellOverflowInto(1, 1), cellOverflowInto(1, 2));
    [numrows,numcols] = size(dem);
    [second_pit, ~] = find(pitId == second_pit_ID);

    spilloverElevation(second_pit) = NaN;
    outletCell(second_pit, :) = [0, 0];
    cellOverflowInto(second_pit, :) = [0, 0];

    if isinf(vca(second_pit))
        % re-ID first pit, raise/fill first pit cells
        raisedPointsCount = 0;
        indicesToCheck = int16(zeros(areaCellCount(1), 2));
        indicesToCheck(1, :) = pitCell(1, :);
        j = 2;
        for i = 1 : size(indicesToCheck, 1)
            r = indicesToCheck(i, 1);
            c = indicesToCheck(i, 2);
            if (dem(r,c) <= spilloverElevation(1))
                raisedPointsCount = raisedPointsCount + 1;
                dem(r,c) = spilloverElevation(1);
            else
                pits(r,c) = pitId(second_pit);
            end

            if (dem(r,c) < spilloverElevation(second_pit))
                volume(second_pit) = volume(second_pit) + ((spilloverElevation(second_pit) - dem(r,c))*cellsize*cellsize);
            end

            for x = 1 : 3 
                for y = 1 : 3
                    if flow_direction_parents(r,c,y,x)
                        indicesToCheck(j, :) = [r+y-2, c+x-2];
                        j = j + 1;
                    end
                end
            end
        end
        
        % Resolve the flow direction of the filled cells.  Begin by directing
        % the spillover point to the pre-identified outletSpillover point
        indicesToCheck = int16(zeros(raisedPointsCount, 2));
        flow_direction(outletCell(1,1), outletCell(1,2), :) = cellOverflowInto(1,:);
        rr = outletCell(1, 1);
        cc = outletCell(1, 2);
        r = cellOverflowInto(1, 1);
        c = cellOverflowInto(1, 2);
        flow_direction_parents(r,c, rr-r+2, cc-c+2) = true;
        pits(outletCell(1, 1), outletCell(1, 2)) = pitId(second_pit);
        indicesToCheck(1, :) = outletCell(1,:);
        j = 2;
        for i = 1 : size(indicesToCheck, 1)
            r = indicesToCheck(i, 1);
            c = indicesToCheck(i, 2);
            for x = -1 : 1 % loop through neighboring cells
                for y = -1 : 1
                    if x == 0 && y ==0 
                        continue; % skip center cell of 3x3 neighborhood
                    end
                    % Check if the point needs to be resolved, but hasn't
                    % already be resolved. This avoids redundant additions of
                    % points to the "next up" list of points to be resolved.
                    if ((pits(r+y, c+x) == pitId(1)) && (pits(r+y, c+x) ~= pitId(second_pit)))
    %                     if flow_direction(r+y, c+x, :) ~= -1
    %                         flow_direction_parents(flow_direction(r+y, c+x, 1), flow_direction(r+y, c+x, 2), (r+y)-(flow_direction(r+y, c+x, 1))+2, (c+x)-(flow_direction(r+y, c+x, 2))+2) = false;
    %                     end
                        flow_direction(r+y, c+x, :) = [r,c];
    %                     flow_direction_parents(r,c,y+2, x+2) = true;
                        pits(r+y, c+x) = pitId(second_pit);
                        indicesToCheck(j, :) = [r+y, c+x];
                        j = j + 1;
                    end
                end
            end
        end

        % Correct flow direction parents
        for i = 1 : size(indicesToCheck, 1)
            r = indicesToCheck(i, 1);
            c = indicesToCheck(i, 2);
            flow_direction_parents(r,c,:,:) = false;
            for x = -1 : 1 % loop through neighboring cells
                for y = -1 : 1
                    if x == 0 && y == 0 % skip center (target) cell of 3x3 neighborhood
                        continue;
                    end

                    if r+y > numrows || r+y < 1 || c+x > numcols || c+x < 1
                        continue; % skip neighbors outside the matrix range
                    end
                    if flow_direction(r+y, c+x, :) ~= 0
                        if isequal(flow_direction(r+y, c+x, 1), r) && isequal(flow_direction(r+y, c+x, 2), c)
                            flow_direction_parents(r, c, y+2, x+2) = true;
                        end
                    end
                end
            end
        end

        areaCellCount(second_pit) = areaCellCount(second_pit) +  areaCellCount(1);
        filledVolume(second_pit) = volume(1) + filledVolume(second_pit);
        
    else
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        % Begin identifying the spillover location in the first pit.
        indicesToCheck = int16(zeros(areaCellCount(1), 2));
        indicesToCheck(1, :) = pitCell(1, :);
        j = 2;
        for i = 1 : size(indicesToCheck, 1)
            r = indicesToCheck(i, 1);
            c = indicesToCheck(i, 2);
            for x = -1 : 1 % loop through neighboring cells
                for y = -1 : 1
                    if x == 0 && y ==0 
                        continue; % skip center cell of 3x3 neighborhood
                    end
                    if r+y > numrows || r+y < 1 || c+x > numcols || c+x < 1
                        continue; % skip neighbors outside the dem range
                    end
                    if isnan(dem(r+y, c+x))
                        continue;
                    end

                    % First, check for neighbors to check next.
                    if flow_direction_parents(r,c,y+2, x+2)
                        indicesToCheck(j, :) = [r+y, c+x];
                        j = j + 1;
                    end

                    % Next, check if it is on the border.
                    if ((pits(r+y, c+x) ~= pitId(second_pit)) && (pits(r+y, c+x) ~= pitId(1)))
                        % if minimum ridge elevation value is still NaN from
                        % initialization or if the ridge elevations (in and the
                        % neighbor just out of the pit) are BOTH less than the current
                        % minimum ridge elevation in the pit
                        cur_cell_elev = dem(r, c);
                        neighbor_elev = dem(r+y, c+x);
                        if isnan(spilloverElevation(second_pit)) || (cur_cell_elev <= spilloverElevation(second_pit) && neighbor_elev <= spilloverElevation(second_pit))
                            spilloverElevation(second_pit) = max([neighbor_elev, cur_cell_elev]);
                            outletCell(second_pit, :) = [r, c];
                            cellOverflowInto(second_pit, :) = [r+y, c+x];
                        end
                    end
                end
            end
        end

        % Extend the search for the spillover locaiton to the second pit. Only 
        % at the end of this step is the spillover location is known.
        indicesToCheck = int16(zeros(areaCellCount(second_pit), 2));
        indicesToCheck(1, :) = pitCell(second_pit, :);
        j = 2;
        for i = 1 : size(indicesToCheck, 1)
            r = indicesToCheck(i, 1);
            c = indicesToCheck(i, 2);
            for x = -1 : 1 % loop through neighboring cells
                for y = -1 : 1
                    if x == 0 && y ==0 
                        continue; % skip center cell of 3x3 neighborhood
                    end
                    if r+y > numrows || r+y < 1 || c+x > numcols || c+x < 1
                        continue; % skip neighbors outside the dem range
                    end
                    if isnan(dem(r+y, c+x))
                        continue;
                    end

                    % First, check for neighbors to check next.
                    if flow_direction_parents(r,c,y+2, x+2)
                        indicesToCheck(j, :) = [r+y, c+x];
                        j = j + 1;
                    end

                    % Next, check if it is on the border.
                    if ((pits(r+y, c+x) ~= pitId(second_pit)) && (pits(r+y, c+x) ~= pitId(1)))
                        % if minimum ridge elevation value is still NaN from
                        % initialization or if the ridge elevations (in and the
                        % neighbor just out of the pit) are BOTH less than the current
                        % minimum ridge elevation in the pit
                        cur_cell_elev = dem(r, c);
                        neighbor_elev = dem(r+y, c+x);
                        if isnan(spilloverElevation(second_pit)) || (cur_cell_elev <= spilloverElevation(second_pit) && neighbor_elev <= spilloverElevation(second_pit))
                            spilloverElevation(second_pit) = max([neighbor_elev, cur_cell_elev]);
                            outletCell(second_pit, :) = [r, c];
                            cellOverflowInto(second_pit, :) = [r+y, c+x];
                        end
                    end
                end
            end
        end

        % Update the filled volumes of the pit
        filledVolume(second_pit) = volume(1) + filledVolume(second_pit);
        volume(second_pit) = filledVolume(second_pit);

        % re-ID first pit, raise/fill first pit cells
        raisedPointsCount = 0;
        indicesToCheck = int16(zeros(areaCellCount(1), 2));
        indicesToCheck(1, :) = pitCell(1, :);
        j = 2;
        for i = 1 : size(indicesToCheck, 1)
            r = indicesToCheck(i, 1);
            c = indicesToCheck(i, 2);
            if (dem(r,c) <= spilloverElevation(1))
                raisedPointsCount = raisedPointsCount + 1;
                dem(r,c) = spilloverElevation(1);
            else
                pits(r,c) = pitId(second_pit);
            end

            if (dem(r,c) < spilloverElevation(second_pit))
                volume(second_pit) = volume(second_pit) + ((spilloverElevation(second_pit) - dem(r,c))*cellsize*cellsize);
            end

            for x = 1 : 3 
                for y = 1 : 3
                    if flow_direction_parents(r,c,y,x)
                        indicesToCheck(j, :) = [r+y-2, c+x-2];
                        j = j + 1;
                    end
                end
            end
        end

        % Calculate volume through the second pit
        indicesToCheck = int16(zeros(areaCellCount(second_pit), 2));
        indicesToCheck(1, :) = pitCell(second_pit, :);
        j = 2;
        for i = 1 : size(indicesToCheck, 1)
            r = indicesToCheck(i, 1);
            c = indicesToCheck(i, 2);
            if (dem(r,c) < spilloverElevation(second_pit))
                volume(second_pit) = volume(second_pit) + ((spilloverElevation(second_pit) - dem(r,c))*cellsize*cellsize);
            end

            for x = 1 : 3 
                for y = 1 : 3
                    if flow_direction_parents(r,c,y,x)
                        indicesToCheck(j, :) = [r+y-2, c+x-2];
                        j = j + 1;
                    end
                end
            end
        end

        % Resolve the flow direction of the filled cells.  Begin by directing
        % the spillover point to the pre-identified outletSpillover point
        indicesToCheck = int16(zeros(raisedPointsCount, 2));
    %     if flow_direction(outletCell(1,1), outletCell(1,2), :) ~= -1
    %         flow_direction_parents(flow_direction(outletCell(1,1), outletCell(1,2), 1), flow_direction(outletCell(1,1), outletCell(1,2), 2), (outletCell(1,1))-(flow_direction(outletCell(1,1), outletCell(1,1), 1))+2, (outletCell(1,2))-(flow_direction(outletCell(1,1), outletCell(1,2), 2))+2) = false;
    %     end
        flow_direction(outletCell(1,1), outletCell(1,2), :) = cellOverflowInto(1,:);
        rr = outletCell(1, 1);
        cc = outletCell(1, 2);
        r = cellOverflowInto(1, 1);
        c = cellOverflowInto(1, 2);
        flow_direction_parents(r,c, rr-r+2, cc-c+2) = true;
        pits(outletCell(1, 1), outletCell(1, 2)) = pitId(second_pit);
        indicesToCheck(1, :) = outletCell(1,:);
        j = 2;
        for i = 1 : size(indicesToCheck, 1)
            r = indicesToCheck(i, 1);
            c = indicesToCheck(i, 2);
            for x = -1 : 1 % loop through neighboring cells
                for y = -1 : 1
                    if x == 0 && y ==0 
                        continue; % skip center cell of 3x3 neighborhood
                    end
                    % Check if the point needs to be resolved, but hasn't
                    % already be resolved. This avoids redundant additions of
                    % points to the "next up" list of points to be resolved.
                    if ((pits(r+y, c+x) == pitId(1)) && (pits(r+y, c+x) ~= pitId(second_pit)))
                        flow_direction(r+y, c+x, :) = [r,c];
                        pits(r+y, c+x) = pitId(second_pit);
                        indicesToCheck(j, :) = [r+y, c+x];
                        j = j + 1;
                    end
                end
            end
        end

        % Correct flow direction parents
        for i = 1 : size(indicesToCheck, 1)
            r = indicesToCheck(i, 1);
            c = indicesToCheck(i, 2);
            flow_direction_parents(r,c,:,:) = false;
            for x = -1 : 1 % loop through neighboring cells
                for y = -1 : 1
                    if x == 0 && y == 0 % skip center (target) cell of 3x3 neighborhood
                        continue;
                    end

                    if r+y > numrows || r+y < 1 || c+x > numcols || c+x < 1
                        continue; % skip neighbors outside the matrix range
                    end
                    if flow_direction(r+y, c+x, :) ~= 0
                        if isequal(flow_direction(r+y, c+x, 1), r) && isequal(flow_direction(r+y, c+x, 2), c)
                            flow_direction_parents(r, c, y+2, x+2) = true;
                        end
                    end
                end
            end
        end

        areaCellCount(second_pit) = areaCellCount(second_pit) +  areaCellCount(1);
        vca(second_pit) = volume(second_pit)/(areaCellCount(second_pit).*(cellsize^2));

        if vca(second_pit) < 0
            vca(second_pit) = Inf;
        end
    end
    
    vca(1) = NaN;
    [~, order] = sort(vca);
    pitId = pitId(order);
    pitCell = pitCell(order, :);
    areaCellCount = areaCellCount(order);
    spilloverElevation = spilloverElevation(order);
    vca = vca(order);
    volume = volume(order);
    filledVolume = filledVolume(order);
    outletCell = outletCell(order, :);
    cellOverflowInto = cellOverflowInto(order, :);
    
    cedarCreekPitId = pits(CedarCreekOutletPoint(1), CedarCreekOutletPoint(2));
    [cedarCreekPitIndex, ~] = find(pitId == cedarCreekPitId);
    number_of_pits(idx) = sum(~isnan(vca));
    runoff_areas(idx) = areaCellCount(cedarCreekPitIndex).*0.0001;
    mean_depth(idx) = nanmean(vca);
    std_depth(idx) = nanstd(vca);
    mean_area(idx) = nanmean(areaCellCount)*0.0001*cellsize^2;
    std_area(idx) = nanstd(areaCellCount)*0.0001*cellsize^2;
    storage_volume(idx) = ((nonNanCount - areaCellCount(cedarCreekPitIndex))*rainfall_excess(idx)*cellsize^2) + filledVolume(cedarCreekPitIndex);
    runoff_volume(idx) = (areaCellCount(cedarCreekPitIndex)*rainfall_excess(idx)*cellsize^2) - filledVolume(cedarCreekPitIndex);
    
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
    
    idx = idx + 1
end
end


