function[times, dem, flow_direction, pits, rainfall_excess, number_of_pits, depthFlow] = ...
    fillDepressions(fillRainfallExcess, dem, flow_direction, pits, pairs, cellIndexes, pitId, pitCell, areaCellCount, spilloverElevation, vca, volume, filledVolume, cellOverflowInto, R, visualize_merging)
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
% number_of_pits - Number of depressions remaining at each iteration.


%First, find the pit containing the Cedar Creek outlet, and set its
%spillover time to infinity.  This will prevent it from overflowing into
%other depressions. Only other depressions may overflow into it.
cellsize = R.CellExtentInWorldX; % DEM cellsize in meters

max_id = max(pitId);
nonNanCount = nansum(nansum(~isnan(dem)));
disp(strcat('Cells with defined values: ', num2str(nonNanCount)));

% Find the total number of iterations/merges
potential_merges = max_id;
% Preallocate arrays of variables that will be used in analysis
rainfall_excess = zeros(potential_merges, 1); % meters
number_of_pits = zeros(potential_merges, 1); 
times = zeros(potential_merges, 1); % s
finalOrder = nan(potential_merges, 1);

% Get initial values
rainfall_excess(1) = 0;
number_of_pits(1) = sum(~isnan(vca));

pitsColormap = rand(max(max(pits))+1, 3); % random color for every depression, plus one for depression ID 0.
pitsColormap(1,:) = 1; % make pit ID 0 white
RGB = ind2rgb(pits, pitsColormap);
image(RGB);
axis equal;
set(gca,'visible','off');
set(gca,'position',[0 0 1 1], 'units', 'normalized');
drawnow;

depthFlow = nan(size(dem));
idx = 2;
fprintf('Filling depressions...')
[~, first_pit] = min(vca);
second_pit = pits(cellOverflowInto(first_pit));

while (vca(first_pit) <= (fillRainfallExcess/1000)) && (idx <= potential_merges+1)
    try
    newTic = tic();
    finalOrder(idx-1) = first_pit;
    rainfall_excess(idx) = vca(first_pit);
    
    % re-ID first pit, raise/fill first pit cells
    lessThanSpillover = cellIndexes{first_pit}(dem(cellIndexes{first_pit}) <= spilloverElevation(first_pit));
    depthFlow(lessThanSpillover) = vca(first_pit); % mark depthFlow with vca
    dem(lessThanSpillover) = spilloverElevation(first_pit); % raise DEM

    %update flow directions
    flow_direction(pitCell(first_pit)) = cellOverflowInto(first_pit);

    pits([cellIndexes{first_pit}]) = second_pit; % update pit matrix
    
    idx
    
    % handle pits now connected to the edge of the DEM
    if (second_pit > 0)
        
        areaCellCount(second_pit) = areaCellCount(first_pit) + areaCellCount(second_pit);

        % Find spillover elevation and location; append pairs, prune, and find
        pairs{second_pit} = [pairs{first_pit}; pairs{second_pit}];
        pairs{second_pit} = pairs{second_pit}(any(pits(pairs{second_pit}) ~= second_pit, 2), :);
        [val, ord] = min(max(dem(pairs{second_pit}(:, 1:2)), [], 2));

        if idx == potential_merges
            times(idx-1) = toc(newTic);
            break;
        end

        spilloverElevation(second_pit) = val;
        cellOverflowInto(second_pit) = pairs{second_pit}(ord, pits(pairs{second_pit}(ord, :)) ~= second_pit);

        % compute filled volume and initialize volume calculation (more below)
        filledVolume(second_pit) = volume(first_pit) + filledVolume(second_pit);
        volume(second_pit) = filledVolume(second_pit);

        cellIndexes{second_pit} = [cellIndexes{second_pit}; cellIndexes{first_pit}];
        lessThans = cellIndexes{second_pit}(dem(cellIndexes{second_pit}) <= spilloverElevation(second_pit));

        volume(second_pit) = volume(second_pit) + sum((spilloverElevation(second_pit) - dem(lessThans) )*cellsize*cellsize);
        lessThans = [];

        vca(second_pit) = volume(second_pit)/((cellsize^2).*areaCellCount(second_pit));
        if vca(second_pit) < 0
            vca(second_pit) = Inf;
        end
    end
  
    vca(first_pit) = NaN;     
    cellIndexes{first_pit} = [];
    pairs{first_pit} = [];
    
    if (visualize_merging) 
        number_of_pits(idx) = sum(~isnan(vca));
        RGB = ind2rgb(pits, pitsColormap);
        image(RGB);
        axis equal;
        set(gca,'visible','off');
        set(gca,'position',[0 0 1 1], 'units', 'normalized');
        drawnow;
    end
    
    idx = idx + 1;
    [~, first_pit] = min(vca);
    second_pit = pits(cellOverflowInto(first_pit));
    times(idx-1) = toc(newTic);
    catch ME
       idx
       break;
    end
end
number_of_pits(idx) = sum(~isnan(vca));
RGB = ind2rgb(pits, pitsColormap);
image(RGB);
axis equal;
set(gca,'visible','off');
set(gca,'position',[0 0 1 1], 'units', 'normalized');
drawnow;
fprintf('Done')
fprintf('\n')

end


