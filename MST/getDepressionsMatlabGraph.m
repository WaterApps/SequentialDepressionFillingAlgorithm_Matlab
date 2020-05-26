function[pits, pairs, cellIndexes, pitId, pitCell, areaCellCount, spilloverElevation, vca, volume, filledVolume, cellOverflowInto, directedGraph] = getDepressionsMatlabGraph(dem, flow_direction, flow_direction_parents, cellsize)
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
edge_pits = int32(find(flow_direction == -2));
pitId = int32(1 : pit_count); % id of each pit
edgePitId = int32(-1 : -1 : (-1*length(edge_pits)-1));
pitCell = int32(find(flow_direction == -1)); % pit bottom cell index
edgePitCell = int32(find(flow_direction == -2));
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
cellIndexes = cell(length(pitId), 1);
for p = 1 : length(pitCell)
    j = 1; % current length of cellIndexes{p}
    i = 1; % iterator
    chunk = 50;
    cellIndexes{p} = nan(chunk, 1);
    cellIndexes{p}(1) = pitCell(p);
    while i <= j
        parents = flow_direction_parents{cellIndexes{p}(i)};
        pits(parents) = p;
        k = j + length(parents);
        if (k > chunk)
            cellIndexes{p}(chunk+1:chunk+50) = zeros(50, 1);
            chunk = chunk + 50;
        end
        cellIndexes{p}(j+1 : k) = parents;
        j = k;
        i = i+1;
    end
    cellIndexes{p}(j+1:end) = [];
    areaCellCount(p) =j;
end

pits(edgePitCell) = length(pitCell) + 1 : length(pitCell) + length(edgePitCell);
for p = 1 : length(edgePitCell)
    j = 1; % current length of cellIndexes{p}
    i = 1; % iterator
    chunk = 50;
    
    indexes = nan(chunk, 1);
    indexes(1) = edgePitCell(p);
    while i <= j
        parents = flow_direction_parents{indexes(i)};
        pits(parents) = length(pitCell) + p;
        k = j + length(parents);
        if (k > chunk)
            indexes(chunk+1:chunk+50) = zeros(50, 1);
            chunk = chunk + 50;
        end
        indexes(j+1 : k) = parents;
        j = k;
        i = i+1;
    end
end
fprintf('\n');
fprintf('Gathering data for each pit...');

indexes = reshape(1: numel(flow_direction), size(flow_direction));
pairs = cell(pit_count, 1);
[numrows, numcols] = size(dem);

adjacents = cell(pit_count, 1);
maxPit = length(edgePitCell) + length(pitCell) + 1;
boundEdges = [pits(edgePitCell), zeros(size(edgePitCell))+maxPit, inf(size(edgePitCell)).*-1, zeros(size(edgePitCell)),zeros(size(edgePitCell))];

for p = 1 : pit_count
    disp(p)
    row = 1;
    l = 1;
    pairs{p} = zeros(length(cellIndexes{p})*8, 2);
    adjacents{p} = [];
    for i = 1 : length(cellIndexes{p}) % walk through indicesToCheck
        [r, c] = ind2sub(size(pits), cellIndexes{p}(i));
        % First, loop over 3x3 neighborhood to find indicesToCheck
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
                neighbor = indexes(r+y, c+x);
                curCell = indexes(r, c);
                % First, check for neighbors to check next.
                if pits(neighbor) ~= pits(curCell)
                    pairs{p}(l, 1) = min([neighbor, curCell]); % pairs are always oriented starting at lowest linear index
                    pairs{p}(l, 2) = max([neighbor, curCell]);
                    
                    % Now handle 
                    maxOverflow = max(dem(neighbor), dem(curCell));
                    % get or create adjacency entry for current pit
                    entry = [];
                    if (isempty(adjacents{p})) 
                        entry = [];
                    else
                        entry = find(adjacents{p}(:, 2) == pits(neighbor));
                    end

                    if (isempty(entry))
                        entry = row;
                        adjacents{p}(row, :) = [pits(curCell), pits(neighbor), curCell, neighbor, inf];
                        row = row + 1;
                    end
                    
                    curMin = adjacents{p}(entry, 5);
                    adjacents{p}(entry, 5) = min([maxOverflow, curMin]);
                    l = l+1;
                end
            end
        end
    end
    
    pairs{p} = pairs{p}(any(pairs{p}, 2), :); % trim zero rows
    pairs{p} = sortrows(pairs{p}, [1,2]);
    [val, ord] = min(max(dem(pairs{p}(:, 1:2)), [], 2)); % get minimum spillover elevation
    spilloverElevation(p) = val;
    cellOverflowInto(p) = pairs{p}(ord, pits(pairs{p}(ord, :)) ~= p);
    
    % Accumulate those elevations that are less than spillover
    lessThans = cellIndexes{p}(dem(cellIndexes{p}) <= spilloverElevation(p));
    volume(p) = sum((spilloverElevation(p) - dem(lessThans))*cellsize*cellsize); 
    lessThans = []; 
    
    filledVolume(p) = 0;

    vca(p) = volume(p)/((cellsize^2).*areaCellCount(p));
    if vca(p) < 0
        vca(p) = Inf;
    end
end


adjacents = vertcat(adjacents{:});
adjacents = [adjacents; boundEdges];

tab = table([adjacents(:, 1), adjacents(:, 2)], adjacents(:,5), adjacents(:, 3), adjacents(:, 4), 'VariableNames', {'EndNodes', 'Weight', 'fromCellIndex', 'toCellIndex'});
G = simplify(graph(tab));
T = minspantree(G);
% Making it directed isn't quite right at the moment
% DG = digraph(T.Edges, T.Nodes);

figure(1)
p = plot(G);
highlight(p, T, 'EdgeColor', 'r', 'NodeColor', 'r')

figure(2)
imagesc(pits);

figure(3);
imagesc(dem);

figure(4);
directedGraph = orientGraph(T);
plot(directedGraph, 'ArrowSize', 8)

disp('done')

end

function[directedGraph] = orientGraph(gr)
pitCount = height(gr.Nodes) -1;
sp = adjacency(gr);

from = zeros(pitCount, 1);
to = zeros(pitCount, 1);
Weights = zeros(pitCount, 1);
queueIdx = 1;
queue = zeros(pitCount, 1);
queue(queueIdx) = pitCount + 1;
fromCellIndex = zeros(pitCount, 1);
toCellIndex = zeros(pitCount, 1);
previous = zeros(pitCount, 1);

previous = 0; % first iteration of 'order' should have no previous
for i = 1 : pitCount
    disp(i)
    order = graphtraverse(sp, queue(i), 'Directed', false, 'Depth', 1);
    if (length(order) < 2) 
        continue;
    end
    for j = 2 : length(order)
        if (order(j) == previous(i))
            continue
        end
        
        to(queueIdx) = order(1);
        from(queueIdx) = order(j);
        edgeIdx = findedge(gr, order(1), order(j));
        Weights(queueIdx) = gr.Edges.Weight(edgeIdx);
        fromCellIndex(queueIdx) = gr.Edges.fromCellIndex(edgeIdx);
        toCellIndex(queueIdx) = gr.Edges.toCellIndex(edgeIdx);
        previous(queueIdx + 1) = order(1);
        queue(queueIdx + 1) = order(j);
        queueIdx = queueIdx + 1;
    end
end
EndNodes = [from, to];
EdgeTable = table(EndNodes, Weights, toCellIndex, fromCellIndex);
directedGraph = digraph(EdgeTable);
end

