function[dem, flow_direction, pits, depthFlow, rainfall_excess, runoff] = ...
    fillDepressionsMovieFlowAccumulation(fillRainfallExcess, dem, flow_direction, pits, pairs, cellIndexes, pitId, pitCell, edgePitCell, areaCellCount, spilloverElevation, vca, volume, filledVolume, cellOverflowInto, perimeterCell, R, visualize_merging)
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
profile on;
cellsize = R.CellExtentInWorldX; % DEM cellsize in meters

max_id = max(pitId);
nonNanCount = nansum(nansum(~isnan(dem)));
disp(strcat('Cells with defined values: ', num2str(nonNanCount)));

% Find the total number of iterations/merges
potential_merges = max_id;
% Preallocate arrays of variables that will be used in analysis
rainfall_excess = zeros(potential_merges, 1); % meters
times = zeros(potential_merges, 1); % s
finalOrder = nan(potential_merges, 1);
runoff = zeros(potential_merges, 1); % area in terms of cells!

% Get initial values
rainfall_excess(1) = 0;
runoff(1) = sum(sum(pits < 0));

ii = double(flow_direction(~isnan(flow_direction) & flow_direction > 0));
jj = find(~isnan(flow_direction) & flow_direction > 0);
flow_direction_parents = cell(size(flow_direction));
for i = 1 : numel(ii)
    flow_direction_parents{ii(i)} = [flow_direction_parents{ii(i)}, jj(i)];
end
flow_accumulation = nan(size(flow_direction));
flow_accumulation(~isnan(flow_direction)) = 0;
for i = 1 : length(ii)
    if flow_accumulation(ii(i)) == 0 % begin recursion if it hasn't been yet visited
        parentSum = recursionThree(ii(i));
    end
end

function tot = recursionThree(i)
    flow_accumulation(i) =  flow_accumulation(i) + 1; % count yourself
    parents = flow_direction_parents{i};
    for parent = parents
        if (flow_accumulation(parent) == 0)
            parentSum = recursionThree(parent);
        else
            parentSum = flow_accumulation(parent);
        end
        flow_accumulation(i) =  flow_accumulation(i) + parentSum;
    end
    tot = flow_accumulation(i);
    return
end
 
%Setup Movie files
vid = VideoWriter('catchments.mp4');
open(vid);
a = figure(6);
n = 1;
imagesc(flow_accumulation, 'AlphaData',~isnan(flow_accumulation));
title('Flow Accumulation');
cmap = colormap('gray'); % base it on grayscale
colormap(cmap);
h = colorbar;
set(gca,'colorscale','log');
axis equal;
axis tight manual;
set(gca,'visible','off');
set(gca,'position',[0 0 1 1], 'units', 'normalized');
drawnow;
F(n) = getframe(a);
writeVideo(vid, F(n))

depthFlow = nan(size(dem));
idx = 2;
fprintf('Filling depressions...')
[~, first_pit] = min(vca);
second_pit = pits(cellOverflowInto(first_pit));

while (vca(first_pit) <= (fillRainfallExcess/1000)) && (idx <= potential_merges+1)
    newTic = tic();
    finalOrder(idx-1) = first_pit;
    rainfall_excess(idx) = vca(first_pit);
    runoff(idx) = sum(sum(pits < 0));
    
    % re-ID first pit, raise/fill first pit cells
    lessThanSpillover = cellIndexes{first_pit}(dem(cellIndexes{first_pit}) <= spilloverElevation(first_pit));
    depthFlow(lessThanSpillover) = vca(first_pit); % mark depthFlow with vca
    dem(lessThanSpillover) = spilloverElevation(first_pit); % raise DEM

    %update flow directions
    %flow_direction(pitCell(first_pit)) = cellOverflowInto(first_pit);
    previousCell = cellOverflowInto(first_pit);
    curCell = perimeterCell(first_pit);
    while flow_direction(curCell) > 0
        nextCell = flow_direction(curCell); % store the old flow direction;
        if (~isempty(flow_direction_parents{curCell}(flow_direction_parents{curCell} == previousCell)))
           flow_direction_parents{curCell}(flow_direction_parents{curCell} == previousCell) = [];
        end
        flow_direction(curCell) = previousCell; % set the new flow direction;
%         if isempty(find(flow_direction_parents{previousCell} == curCell, 1))
            flow_direction_parents{previousCell} = [flow_direction_parents{previousCell}, curCell];
%         end
        previousCell = curCell; % advance for the next iteration
        curCell = nextCell; % advance for the next iteration
    end
    if first_pit == 108
       curCell
    end
    flow_direction(curCell) = previousCell;
    if (~isempty(flow_direction_parents{curCell}(flow_direction_parents{curCell} == previousCell)))
       flow_direction_parents{curCell}(flow_direction_parents{curCell} == previousCell) = [];
    end
    flow_direction_parents{previousCell} = [flow_direction_parents{previousCell}, curCell];
    
    
    pits([cellIndexes{first_pit}]) = second_pit; % update pit matrix

    flow_accumulation(pits == second_pit) = 0;
    
    % handle pits not connected to the edge of the DEM
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
        
        %Now fix flow accumulation for the pits involved in the spillover
        parentSum = recursionThree(pitCell(second_pit));
        if (sum(flow_accumulation([cellIndexes{first_pit}]) == 0) > 0)
            sum(flow_accumulation([cellIndexes{first_pit}]) == 0)
        elseif (sum(flow_accumulation([cellIndexes{second_pit}]) == 0) > 0)
            sum(flow_accumulation([cellIndexes{second_pit}]) == 0)
        end
        
    else
        %Now fix flow accumulation for the pits involved in the spillove\
        parentSum = recursionThree(edgePitCell(-1*second_pit));
        if (sum(flow_accumulation([cellIndexes{first_pit}]) == 0) > 0)
            sum(flow_accumulation([cellIndexes{first_pit}]) == 0)
        elseif (sum(flow_accumulation(pits == second_pit) == 0) > 0)
            sum(flow_accumulation(pits == second_pit) == 0)
        end
    end
    
    vca(first_pit) = NaN;
    cellIndexes{first_pit} = [];
    pairs{first_pit} = [];
    
    if (visualize_merging) && rainfall_excess(idx) > 0.018
        imagesc(flow_accumulation, 'AlphaData',~isnan(flow_accumulation));
        title('Flow Accumulation');
        cmap = colormap('gray'); % base it on grayscale
        colormap(cmap);
        h = colorbar;
        set(gca,'colorscale','log');
        axis equal;
        axis tight manual;
        set(gca,'visible','off');
        set(gca,'position',[0 0 1 1], 'units', 'normalized');
        drawnow;
        F(n) = getframe(a);
        writeVideo(vid, F(n))
        n = n + 1;
    end
    
    idx = idx + 1
    [~, first_pit] = min(vca);
    second_pit = pits(cellOverflowInto(first_pit));
    times(idx-1) = toc(newTic);

end
% RGB = ind2rgb(pits, pitsColormap);
% image(RGB);
% axis equal;
% set(gca,'visible','off');
% set(gca,'position',[0 0 1 1], 'units', 'normalized');
drawnow;
fprintf('Done')
fprintf('\n')

end


