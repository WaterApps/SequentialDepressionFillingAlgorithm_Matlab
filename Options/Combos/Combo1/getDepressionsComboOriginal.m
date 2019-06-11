function[pits, pitId, pitCell, areaCellCount, spilloverElevation, vca, volume, filledVolume, cellOverflowInto] = getDepressionsComboOriginal(dem, flow_direction, cellsize)
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

% Fill out the pit ID image/matrix
pits(~isnan(dem)) = 0;
pits(pitCell) = pitId; % Mark the pit bottom in the pits image. This is 
% used in recursiveMarkPits to tell what ID to mark each cell during tail
% recursion.
for p = 1 : numel(pits)
    if (isnan(pits(p)))
        continue
    end
    if (pits(p) == 0)
        recursiveMarkPits(p);
    end
    curPitId = pits(p);
    if ~isnan(curPitId)
        areaCellCount(curPitId) = areaCellCount(curPitId) + 1;
    end
end

function childId = recursiveMarkPits(i)
    child = flow_direction(i);
    if pits(child) == 0 % downstream is yet unidentified in the pits matrix
        childId = recursiveMarkPits(child);
        pits(i) = childId;
    else % downstream already set
        childId = pits(child);
        pits(i) = childId;
    end
    return
end


fprintf('\n');
fprintf('Gathering data for each pit...');

[numrows, numcols] = size(pits);
indexes = reshape(1: numel(flow_direction), size(flow_direction));
for p = 1 : length(pitId)
    spilloverElevation(p) = realmax('double');
    indicesToCheck = int32(zeros(areaCellCount(p), 2));
    [r, c] = ind2sub(size(dem), pitCell(p));
    indicesToCheck(1, :) = [r, c];
    j = 2; % index for appending new indicesToCheck
    k = 1;
    for i = 1 : areaCellCount(p)
        r = indicesToCheck(i, 1);
        c = indicesToCheck(i, 2);
        % First, loop over 3x3 neighborhood to find indicesToCheck
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
                % First, check for neighbors to check next.
                if flow_direction(neighbor) == curCell% using the sparse here should be both fast on lookup and low on memory usage
                    indicesToCheck(j, :) = [r+y, c+x];
                    j = j + 1;
                end
                if pits(neighbor) ~= pits(curCell)
                    if isnan(dem(neighbor))
                        continue;
                    end
                    elev = max(dem([neighbor, curCell])); % lines are always oriented starting at lowest linear index
                    if (elev < spilloverElevation(p))
                        spilloverElevation(p) = elev;
                        cellOverflowInto(p) = neighbor;
                    end
                end
            end
        end
    end

    % Compute Volume
    filledVolume(p) = 0;
    indicesToCheck = zeros(areaCellCount(p), 2);
    [r, c] = ind2sub(size(flow_direction), pitCell(p));
    indicesToCheck(1, :) = [r, c];
    j = 2;
    i = 1;
    while i < j
        r = indicesToCheck(i, 1);
        c = indicesToCheck(i, 2);
        i = i + 1;
        volume(p) = volume(p) + (spilloverElevation(p) - dem(r,c))*cellsize*cellsize;
        for x = -1 : 1 % loop through neighboring cells
            for y = -1 : 1
                if x == 0 && y ==0
                    continue; % skip center cell of 3x3 neighborhood
                end

                if r+y > numrows || r+y < 1 || c+x > numcols || c+x < 1
                    continue; % skip neighbors outside the dem range
                end

                curCell = indexes(r, c);
                if flow_direction(r+y, c+x) == curCell && dem(r+y, c+x) <= spilloverElevation(p)
                    indicesToCheck(j, :) = [r+y, c+x];
                    j = j + 1;
                end
            end
        end
    end
    
    vca(p) = volume(p)/((cellsize^2).*areaCellCount(p));
    if vca(p) < 0
        vca(p) = Inf;
    end
end
end