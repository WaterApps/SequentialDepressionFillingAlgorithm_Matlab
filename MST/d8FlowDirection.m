function [flow_direction] = d8FlowDirection(dem)
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
weights = [sqrt(2), 1, sqrt(2); 1 0 1; sqrt(2), 1 sqrt(2)];
flow_direction = nan(size(dem, 1), size(dem, 2));

indexes = reshape(1: numel(flow_direction), size(flow_direction));
[numrows, numcols] = size(dem);
fprintf('\nProcessing flow direction... ')

parfor i = 1 : numel(flow_direction)
    if isnan(dem(i))
        continue
    end
    [r,c] = ind2sub(size(dem), i);
    maxr = r + 1;
    minr = r - 1;
    maxc = c + 1;
    minc = c - 1;
    wminr = 1;
    wmaxr = 3;
    wminc = 1;
    wmaxc = 3;
    
    % Handle border cells 
    if r == numrows
        maxr = numrows;
        wmaxr = 2;
    elseif r == 1
        minr = 1;
        wminr = 2;
    end

    if c == numcols
        maxc = numcols;
        wmaxc = 2;
    elseif c == 1
        minc = 1;
        wminc = 2;
    end
    
    a = ones(maxr-minr+1, maxc-minc+1).*dem(r,c);
    b = dem(minr:maxr, minc:maxc);
    dwd =(a-b)./weights(wminr:wmaxr, wminc:wmaxc);
    dwdMax = max(dwd(:));
    if dwdMax <= 0        
        % TODO: mark edge pits, they shouldn't overflow back into the mix
        % of things, we don't know their true outcome.
        if r == numrows || r == 1 || c == numcols || c == 1
            % mark edge pis as infinite sinks 
            flow_direction(i) = -2;
            continue;
        else
            flow_direction(i) = -1;
            continue;
        end
    end            
    indices = indexes(minr:maxr, minc:maxc);
    found = find(dwd == dwdMax);

    if length(found) == 1
        flow_direction(i) = indices(found);
    else
        % multiple neighboring cells with equivalent downslopes; choose one at random
        flow_direction(i) = indices(datasample(found, 1));
    end
end
fprintf('Done')
fprintf('\n')
flow_direction(isnan(dem)) = NaN;
end
