function[flow_direction] = updateFlowDirectionAdjacent(areaCellCount, pitId, cellOverflowInto, spilloverElevation, flow_direction, dem, pits, indexes)
    % reorient filled zone back toward spillover cell, change pit ID to
    % second_pit, and raise DEM
    [numrows, numcols] = size(flow_direction);
    visited = false(numrows, numcols);
    indicesToCheck = zeros(areaCellCount+1, 2);
    [r, c] = ind2sub(size(flow_direction), cellOverflowInto);
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
                if pits(neighbor) == pitId && dem(neighbor) <= spilloverElevation && ~visited(neighbor)
                    visited(neighbor) = true;
                    j = j + 1;
                    indicesToCheck(j, :) = [r+y, c+x];
                    % add new flow direction parents
                    flow_direction(neighbor) = curCell;
%                     if (neighbor == pitCell)
%                         return;
%                     end
                end
            end
        end
    end
end