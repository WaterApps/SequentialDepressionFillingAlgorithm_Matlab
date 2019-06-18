function[spilloverElevation, cellOverflowInto] = findSpilloverFlowDirection(areaCellCount, pitCell, dem, flow_direction, pits, indexes)
    spilloverElevation = realmax('double');
    indicesToCheck = zeros(areaCellCount, 2);
    [r, c] = ind2sub(size(dem), pitCell);
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
                    if elev < spilloverElevation
                        spilloverElevation = elev;
                        cellOverflowInto = neighbor;
                    end
                end
            end
        end
    end
end


