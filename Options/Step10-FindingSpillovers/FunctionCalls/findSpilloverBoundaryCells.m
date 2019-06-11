function[spilloverElevation, cellOverflowInto, boundaryCells] = findSpilloverBoundaryCells(dem, pits, boundaryCells, indexes)
    % find the spillover of the new combined depression; prune out those
    % which are no longer on the boundary
    spilloverElevation = realmax('double');
    [numrows, numcols] = size(dem);
    for p = length(boundaryCells) : -1 : 1 % reverse order so that we can prune in the same pass
        stillBoundary = false;
        [r, c] = ind2sub(size(dem), boundaryCells(p));
        for x = -1 : 1 % loop through neighboring cells
            for y = -1 : 1
                if x == 0 && y ==0
                    continue; % skip center cell of 3x3 neighborhood
                end
                if r+y > numrows || r+y < 1 || c+x > numcols || c+x < 1
                    continue; % skip neighbors outside the dem range
                end

                neighbor = indexes(r+y, c+x);
                curCell = boundaryCells(p);

                % First, check for neighbors to check next.
                if pits(neighbor) ~= pits(curCell)
                    stillBoundary = true;
                    if isnan(pits(neighbor))
                        continue;
                    end
                    elev = max(dem([neighbor, curCell]));
                    if (elev < spilloverElevation)
                        spilloverElevation = elev;
                        cellOverflowInto = neighbor;
                    end
                end
            end
        end
        if (~stillBoundary)%prune out boundaries no longer on the perimeter
            boundaryCells(p) = [];
        end
    end
end
