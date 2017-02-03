function[pits, pitId, pitCell, areaCellCount, spilloverElevation, vca, volume, filledVolume, outletCell, cellOverflowInto] = getDepressions(dem, flow_direction, flow_direction_parents, cellsize)
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
% outletCell - The row and column of the cell where water will begin to
% overflow once the depression is filled.
%
% cellOverflowInto - The row and column of the cell in an adjacent
% depression where water will begin to overflow into.

[numrows, numcols] = size(dem);

% Preallocate variables
pits = int32(zeros(size(dem), 'single')); % matrix identifying pits (pit map)
pit_count = sum(nansum(flow_direction(:,:,1) < 0));
pitId = int32(zeros(pit_count, 1));
pitCell = int16(zeros(pit_count, 2));
areaCellCount = zeros(pit_count, 1); % integer number of cells
spilloverElevation = zeros(pit_count, 1); % meters
vca = zeros(pit_count, 1); % hours
volume = zeros(pit_count, 1); % cubic meters
filledVolume = zeros(pit_count, 1); % cubic meters
outletCell = int16(zeros(pit_count, 2));
cellOverflowInto = int16(zeros(pit_count, 2));

% indicesDrainingToPit(1:pit_count) = struct('allIndices', []);
max_pit_id = 1; % initialize pit numbering
% min_pit_id = - 1;

% Pits must first be identified in the pit matrix in order to return the
% correct pit ID that each pit flows into (if not, many of these pits will
% flow into yet unidentified pits that have ID 0).
pit_counter = 1;
for c = 1 : size(pits, 2)
    for r = 1 : size(pits, 1)
        % Identify positive depressions that require filling
        if flow_direction(r, c, :) < 0
            % Add current cell to list of pit indices
            pitId(pit_counter) = max_pit_id;
            pitCell(pit_counter, :) = [r, c];
            % Identify those cells which flow into the pit
            %preallocate 25
            imax = 25;
            cells_draining_to_pit = int16(zeros(imax, 2));
            %set the first cell to check
            cells_draining_to_pit(1, :) = [r, c];
            i = 1;
            j = 2;
            while cells_draining_to_pit(i, :) ~= 0
                rr = cells_draining_to_pit(i, 1);
                cc = cells_draining_to_pit(i, 2);
                pits(rr, cc) = max_pit_id;
                i = i + 1;
                for x = 1 : 3
                    for y = 1 : 3
                        if flow_direction_parents(rr, cc, y, x)
                            cells_draining_to_pit(j, :) = [rr+y-2, cc+x-2];
                            j = j + 1;
                            if j > imax
                                imax = imax + 25;
                                cells_draining_to_pit(imax, :) = 0;
                            end
                        end
                    end
                end
            end
            areaCellCount(pit_counter) = i-1;
            max_pit_id = max_pit_id + 1;
            pit_counter = pit_counter + 1;

        % Identify negative depressions that are connected to DEM boundary    
%         elseif flow_direction(r,c,:) == 0
%             % Identify those cells which flow into the negative pit, mark them,
%             % but don't retain any data
%             imax = 25;
%             cells_draining_to_pit = int16(zeros(imax, 2));
%             %set the first cell to check
%             cells_draining_to_pit(1, :) = [r, c];
%             i = 1;
%             j = 2;
%             while cells_draining_to_pit(i, :) ~= 0
%                 rr = cells_draining_to_pit(i, 1);
%                 cc = cells_draining_to_pit(i, 2);
%                 pits(rr, cc) = min_pit_id;
%                 i = i + 1;
%                 for x = 1 : 3
%                     for y = 1 : 3
%                         if flow_direction_parents(rr, cc, y, x)
%                             cells_draining_to_pit(j, :) = [rr+y-2, cc+x-2];
%                             j = j + 1;
%                             if j > imax
%                                 imax = imax + 25;
%                                 cells_draining_to_pit(imax, :) = 0;
%                             end
%                         end
%                     end
%                 end
%             end
%             min_pit_id = min_pit_id - 1;
        elseif isnan(flow_direction(r,c, :))
            pits(rr, cc) = 0;
        end
    end
end

% Gather pit_data for each pit
for p = 1 : pit_count
    if pitId(p) > 0
        spilloverElevation(p) = NaN;
        outletCell(p, :) = NaN;
        cellOverflowInto(p, :) = NaN;
        % Find first pit borders for the purpose of 
        indicesToCheck = int16(zeros(areaCellCount(p), 2));
        indicesToCheck(1, :) = pitCell(p, :);
        j = 2;
        for i = 1 : size(indicesToCheck, 1)
            r = indicesToCheck(i, 1);
            c = indicesToCheck(i, 2);
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
                    
                    % First, check for neighbors to check next.
                    if flow_direction_parents(r,c,y+2, x+2)
                        indicesToCheck(j, :) = [r+y, c+x];
                        j = j + 1;
                    end
            
                    % If the neighbor is outside the pit, the border has been
                    % reached.
                    if pits(r+y, c+x) ~= pits(r, c)
                        % if minimum ridge elevation value is still NaN from
                        % initialization or if the ridge elevations (in and the
                        % neighbor just out of the pit) are BOTH less than the current
                        % minimum ridge elevation in the pit
                        cur_cell_elev = dem(r, c);
                        neighbor_elev = dem(r+y, c+x);
                        if isnan(spilloverElevation(p)) || (cur_cell_elev <= spilloverElevation(p) && neighbor_elev <= spilloverElevation(p))
                            spilloverElevation(p) = max([neighbor_elev, cur_cell_elev]);
                            outletCell(p, :) = [r,c];
                            cellOverflowInto(p, :) = [r+y, c+x];
                        end
                    end
                end
            end
        end

        %% Calculate Volume
        volume(p) = 0;
        indicesToCheck = int16(zeros(areaCellCount(p), 2));
        indicesToCheck(1, :) = pitCell(p, :);
        j = 2;
        for i = 1 : size(indicesToCheck, 1)
            r = indicesToCheck(i, 1);
            c = indicesToCheck(i, 2);
            if (dem(r, c) < spilloverElevation(p))
            	volume(p) = volume(p) + ((spilloverElevation(p) - dem(r,c))*cellsize*cellsize);    
            end
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
                    % First, check for neighbors to check next.
                    if flow_direction_parents(r, c, y+2, x+2)
                        indicesToCheck(j, :) = [r+y, c+x];
                        j = j + 1;
                    end
                end
            end
        end

        filledVolume(p) = 0;     
        vca(p) = volume(p)/((cellsize^2).*areaCellCount(p));
        if vca(p) < 0
            vca(p) = Inf;
        end
    
        
    elseif pit_data(p).pitId < 0
        
    end
end
end
