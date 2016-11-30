function [flow_direction, flow_direction_parents] = d8FlowDirectionDrainage(dem, drainage, intensity)
% This function evaluates the DEM and determines the D8 flow direction for
% each cell. Instead of using the powers of two for each of the eight
% neighors or simply 1-8, an angle was used such that it may be possible in
% the future to implement D-Infinity or other multiple outflow direction
% models without as much effort.

% Initialize flow direction as a matrix of -3s. Borders will become NaNs,
% pits resulting from flow direction will be -1, pits resulting from
% excessive drainage rates will be denoted -2 in case this exception
% becomes important, and valid flow directions will range from 0 to 2pi
% (initializing a matrix full of zeros or NaNs may make tests as to whether
% a cell is part of a pit or has yet to be visited unclear if these values
% have multiple identities).

flow_direction = int16(zeros(size(dem, 1), size(dem, 2), 2));
flow_direction_parents = false(size(dem,1), size(dem,2), 3, 3);
[numrows, numcols] = size(dem);
for r = 1 : size(flow_direction, 1)
    for c = 1 : size(flow_direction, 2)
        % Set border cells to NaNs
        if r == numrows || r == 1 || c == numcols || c == 1
            flow_direction(r, c, :) = 0;
            continue;
        end
        % If the cell is draining faster than accumulation then the cell is
        % a pit (given a -2 to denote this special case).
        if drainage(r, c) >= intensity
            flow_direction(r, c, :) = -2;
            continue;
        end
        % If the cell is otherwise NaN, assign a NaN flow direction and
        % continue
        if isnan(dem(r, c))
            flow_direction(r, c, :) = NaN;
            continue;
        end
        max_dist_weighted_drop = NaN;
        % Verify the current_cell hasn't previously been set. May've been set
        % in border or drainage conditions, or by resolveFlatD8FlowDirection if
        % a flat area is located.
        if flow_direction(r, c, :) == 0
            % Compute slope to each neighbor and find minimum
            for x = -1 : 1 % loop through neighboring cells
                for y = -1 : 1
                    if x == 0 && y == 0 % skip center (target) cell of 3x3 neighborhood
                        continue;
                    end
                    % Convert from cartesian to polar coordinates to get flow
                    % direction angle (radians) and distance from center cell to
                    % neighbor (1 for top,bot,left,right, and sqrt(2) for
                    % diagonals).
                    if isnan(dem(r+y, c+x))
                        continue;
                    end
                    distance = pdist([r+y, c+x; r, c], 'euclidean');
                    dist_weighted_drop = (dem(r+y, c+x) - dem(r,c))/distance;

                    % Maintain current minimum slope
                    if (isnan(max_dist_weighted_drop) || dist_weighted_drop < max_dist_weighted_drop) % nan on first iteration
                        max_dist_weighted_drop = dist_weighted_drop;
                        flow_direction(r, c, :) = [r+y, c+x];
                    end
                end
            end

            % Identify flat areas that have one or more neighbors with a
            % distance-weighted drop less than or equal to zero.
            if max_dist_weighted_drop >= 0
                flow_direction(r, c, :) = -1;
            end
        end
    end
end

%Now find parents
for r = 1 : size(flow_direction, 1)
    for c = 1 : size(flow_direction, 2)
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
end

end
