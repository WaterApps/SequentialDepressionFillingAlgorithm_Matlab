function[pits, pitId, pitCell, areaCellCount, spilloverElevation, vca, volume, filledVolume, cellOverflowInto, cellIndexes] = getDepressionsParents(dem, flow_direction, flow_direction_parents, cellsize)
% Generate the matrix identifying depressions in the DEM as well as their
% parameters. 
%
% Inputs: 
%
% dem - elevation matrix
%
% flow_direction - flow direction matrix
%
% flow_direction_parents - matrix of flow direction parents
%
% cellsize - side length of each cell, i.e., DEM resolution (meters).
%
% Outputs (each of the following is an array having a length of the number of
% depressions found in the flow direction matrix operation. Each row
% corresponds to a different depression. Each of these parameters are 
% needed in order to provide consistent mergers of depressions in the 
% filling process):
%
% pitId - The ID number assigned to the depression.
% pitCell - The row and column of the pit cell found in the flow direction 
% matrix.
%
% areaCellCount - Total contributing area expressed as the integer count of
% cells. Multiply by the DEM cell area to get area in square meters.
% spilloverElevation - The elevation which sets the maximum water level of 
% the depression after which it will overflow.
%
% vca - Volume to Contributing Area (VCA) ratios for each depression, also 
% equivalent to the depth of rainfall excess the depression can retain. 
% This value is used to sort the depressions so that they may be filled in 
% order. Computed using volume in cubic meters and area in square meters. 
%
% volume - The retention volumes for each depression (m^3).
%
% filledVolume - The volume (m^3) filled in previous iterations.
%
% cellOverflowInto - The row and column of the cell in an adjacent
% depression where water will begin to overflow into.

% Preallocate variables
disp('Creating Pit Database');

pits = nan(size(dem)); % matrix identifying pits
pit_count = sum(nansum(flow_direction == -1));
pitId = int32(1 : pit_count); % id of each pit
pitCell = int32(find(flow_direction == -1)); % pit bottom cell index
areaCellCount = zeros(pit_count, 1); % integer number of cells
spilloverElevation = zeros(pit_count, 1); % meters
vca = zeros(pit_count, 1); % volume to contributing area ratio (hours)
volume = zeros(pit_count, 1); % (cubic meters)
filledVolume = zeros(pit_count, 1); % (cubic meters)
cellOverflowInto = int32(zeros(pit_count, 1)); %cell index of overflow location

% Pits must first be identified in the pit matrix in order to return the
% correct pit ID that each pit flows into (if not, many of these pits will
% flow into yet unidentified pits that have ID 0).
fprintf('Creating pit ID matrix...')
fprintf('Done');
fprintf('\n');
fprintf(strcat('Number of depressions: ',num2str(length(pitCell))))
fprintf('\n');

% Mark the pit bottom in the pits image
pits(pitCell) = pitId;

cellIndexes = cell(pit_count, 1);
parfor p = 1 : pit_count
    cellIndexes{p} = zeros(100, 1);
    cellIndexes{p} = graphtraverse(flow_direction_parents, pitCell(p), 'Method', 'BFS');
    areaCellCount(p) = nnz(cellIndexes{p});
end
for p = 1 : pit_count
    pits(cellIndexes{p}) = p;
end


parfor p = 1 : length(pitId)
    if pitId(p) > 0
        shape = pits == pitId(p);
        innerPerim = logical(bwperim(shape, 8));
        outerPerim = (conv2(single(innerPerim), ones(3,3), 'same') > 0) & ~shape;
        nanMins = ones(size(shape)).*realmax('double');
        nanMins(innerPerim) = dem(innerPerim);
        %apply a minimum filter in a moving 3x3 window. this applies the
        %minimum value of the 3x3 window to the outerPerimeter. The min
        %between dem values and matlabs realmax is going to resolve to the
        %dem value.
        %Be careful of zeros along the edge of the result of ordfilt2
        nanMins = ordfilt2(nanMins, 1, ones(3));
        nanMins(isnan(dem)) = NaN;
        innerMin = zeros(size(dem));
        innerMin(outerPerim) = nanMins(outerPerim);
        outerMin = zeros(size(dem));
        outerMin(outerPerim) = dem(outerPerim);
        maxCombo = max(innerMin, outerMin);

        % only mark the outer cell;
        spilloverElevation(p) = min(maxCombo(maxCombo > 0));        
        cellOverflowInto(p) = datasample(find(maxCombo==spilloverElevation(p)), 1);
        volume(p) = cellsize*cellsize*sum(spilloverElevation(p) - dem(dem < spilloverElevation(p) & pits == pitId(p)));

        filledVolume(p) = 0;     
        vca(p) = volume(p)/((cellsize^2).*areaCellCount(p));

        if vca(p) < 0
            vca(p) = Inf;
        end
    end
end

% fprintf('\n');
% fprintf('Gathering data for each pit...');
% 
% indexes = reshape(1: numel(flow_direction), size(flow_direction));
% lineSegments = cell(pit_count, 1);
% [numrows, numcols] = size(dem);
% for p = 1 : length(cellIndexes)
%     lineSegments{p} = zeros(length(cellIndexes{p})*8, 3);
%     l = 1;
%     for i = 1 : length(cellIndexes{p}) % walk through indicesToCheck
%         [r, c] = ind2sub(size(pits), cellIndexes{p}(i));
%         % First, loop over 3x3 neighborhood to find indicesToCheck
%         for x = -1 : 1 % loop through neighboring cells
%             for y = -1 : 1
%                 if x == 0 && y ==0
%                     continue; % skip center cell of 3x3 neighborhood
%                 end
%                 if r+y > numrows || r+y < 1 || c+x > numcols || c+x < 1
%                     continue; % skip neighbors outside the dem range
%                 end
%                 if isnan(dem(r+y, c+x))
%                     continue;
%                 end
%                 neighbor = indexes(r+y, c+x);
%                 curCell = indexes(r, c);
%                 % First, check for neighbors to check next.
%                 if pits(neighbor) ~= pits(curCell)
%                     lineSegments{p}(l, 1) = min([neighbor, curCell]); % lines are always oriented starting at lowest linear index
%                     lineSegments{p}(l, 2) = max([neighbor, curCell]);
%                     lineSegments{p}(l, 3) = neighbor; %always spillover to neighboring cell
%                     l = l+1;
%                 end
%             end
%         end
%     end
%     
%     lineSegments{p} = lineSegments{p}(any(lineSegments{p}, 2), :); % trim zero rows
%     lineSegments{p} = sortrows(lineSegments{p}, [1,2]);
%     [val, ord] = min(max(dem(lineSegments{p}(:, 1:2)), [], 2)); % get minimum spillover elevation
%     spilloverElevation(p) = val;
%     cellOverflowInto(p) = lineSegments{p}(ord, 3);
%     
%     % Accumulate those elevations that are less than spillover
%     lessThans = cellIndexes{p}(dem(cellIndexes{p}) <= spilloverElevation(p));
%     volume(p) = sum((spilloverElevation(p) - dem(lessThans))*cellsize*cellsize); 
%     lessThans = []; 
%     
%     filledVolume(p) = 0;
% 
%     vca(p) = volume(p)/((cellsize^2).*areaCellCount(p));
%     if vca(p) < 0
%         vca(p) = Inf;
%     end
% end
end