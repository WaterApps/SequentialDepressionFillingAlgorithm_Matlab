function [flow_direction] = d8FlowDirectionAdjacencyMatrixOperations(dem)
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
newTic = tic();

[r, c] = size(dem);
diagWeight = sqrt(2);
diagVec1 = repmat([ones(c-1, 1); 0], r, 1);% Make the first diagonal vector (for vertical connections given my indexes variable)
sp1 = spdiags(diagVec1, [1], numel(dem), numel(dem));
diagVec2 = [zeros((r*c)-c*(r-1), 1); diagVec1(1:(c*(r-1))).*diagWeight];% Make the second diagonal vector (for diagonal connections)
sp2 = spdiags(diagVec2, [c-1], numel(dem), numel(dem));
diagVec3 = [zeros(c, 1); ones(c*(r-1), 1)]; % Make the third diagonal vector (for horizontal connections)
sp3 = spdiags(diagVec3, c, numel(dem), numel(dem));
diagVec4 = [zeros(2,1); diagVec2(2:end-1)]; % Make the fourth diagonal vector (for anti-diagonal connections)
sp4 = spdiags(diagVec4, c+1, numel(dem), numel(dem));
adj = sp1 + sp2 + sp3 + sp4; % Add the diagonals to a zero matrix
adj = adj + adj.'; % Add the matrix to a transposed copy of itself to make it symmetrics
adj(:, find(isnan(dem))) = 0;
adj(find(isnan(dem)), :) = 0;
% use adjacency matrix to compute adjacency cell array
newadj = cell(size(dem));
nweights = cell(size(dem));
[ii, jj, ss] = find(adj);
for i = 1 : length(ii)
    newadj{ii(i)} = [newadj{ii(i)}, jj(i)];
    nweights{ii(i)} = [nweights{ii(i)}, ss(i)];
end

% weights = [sqrt(2), 1, sqrt(2); 1 0 1; sqrt(2), 1 sqrt(2)];
flow_direction = int64(zeros(size(dem, 1), size(dem, 2)));

% indexes = reshape(1: numel(flow_direction), size(flow_direction));
% [numrows, numcols] = size(dem);
fprintf('\nProcessing flow direction... ')

parfor i = 1 : numel(flow_direction)
    neighbors = newadj{i};
    weights = nweights{i};
    k = length(neighbors);
    if k == 0
       continue; 
    end
    a = dem(i);
    b = dem(neighbors);
    dwd =(a-b)./weights;
    [dwdMax, idx]  = max(dwd(:));
    if dwdMax <= 0
        flow_direction(i) = -1;
        continue;
    end           

    flow_direction(i) = neighbors(idx);
end
toc(newTic)
fprintf('Done')
fprintf('\n')

end

