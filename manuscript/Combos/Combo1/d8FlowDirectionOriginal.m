function [flow_direction] = d8FlowDirectionOriginal(dem)
% Generate the d8 flow direction of water across the DEM surface.
%
% Inputs:
% dem - Matrix of elevations. Has dimensions MxN.
%
% Outputs:
% flow_direction - For each cell, the flow_direction matrix specifies the
% row, column indices of the cell to which water will flow. If the cell is
% a local minima of the 3x3 neighborhood, it will be assigned row and
% column values of -1. Has dimensions MxNx2.
% 
% flow_direction_parents - Boolean matrix specifying whether each of the 
% neighboring cells flow into the current cell. Has size MxNx3x3. For every
% cell (m,n), a 3x3 matrix of neighbors is assigned boolean values based
% on whether or not the neighboring cell flows into the current cell.


% This function evaluates the DEM and determines the D8 flow direction for
% each cell. Instead of using the powers of two for each of the eight
% neighors or simply 1-8, flow directions are encoded as the index of the 
% cell to which flow is directed.

flow_direction = nan(size(dem, 1), size(dem, 2));
[numrows, numcols] = size(dem);
indexes = reshape(1: numel(flow_direction), size(flow_direction));

for r = 1 : size(flow_direction, 1)
    for c = 1 : size(flow_direction, 2)
        % If the cell is NaN in the dem, assign a NaN flow direction and
        % continue
        if isnan(dem(r, c))
            continue;
        end
        max_dist_weighted_drop = realmax('double');

        % Compute slope to each neighbor and find minimum
        for x = -1 : 1 % loop through neighboring cells
            for y = -1 : 1
                if x == 0 && y == 0 % skip center (target) cell of 3x3 neighborhood
                    continue;
                end
                if r+y > numrows || r+y < 1 || c+x > numcols || c+x < 1
                    continue; % skip neighbors outside the dem range
                end
                % Convert from cartesian to polar coordinates to get the
                % distance from center cell to neighbor
                % (1 for top,bot,left,right, and sqrt(2) for diagonals).
                
                distance = pdist([r+y, c+x; r, c], 'euclidean');
                dist_weighted_drop = (dem(r+y, c+x) - dem(r,c))/distance;
                
                % Maintain current minimum slope
                if (dist_weighted_drop < max_dist_weighted_drop) % nan on first iteration
                    max_dist_weighted_drop = dist_weighted_drop;
                    flow_direction(r, c) = indexes(r+y, c+x);
                end
            end

            % Identify flat areas that have one or more neighbors with a
            % distance-weighted drop less than or equal to zero.
            if max_dist_weighted_drop >= 0
                flow_direction(r, c) = -1;
            end
        end
    end
end
fprintf('Done')
fprintf('\n')

end

