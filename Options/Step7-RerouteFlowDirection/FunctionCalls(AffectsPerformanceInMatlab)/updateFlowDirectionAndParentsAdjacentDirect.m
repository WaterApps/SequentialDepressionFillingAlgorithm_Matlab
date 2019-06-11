function[flow_direction, flow_direction_parents] = updateFlowDirectionAndParentsAdjacentDirect(cellOverflowInto, perimeterCell, flow_direction, flow_direction_parents)
    % reorient filled zone back toward spillover cell, change pit ID to
    % second_pit, and raise DEM
    previousCell = cellOverflowInto;
    curCell = perimeterCell;
    while flow_direction(curCell) > 0
        nextCell = flow_direction(curCell); % store the old flow direction;
        
        % remove current cell from array of flow direction parent of the
        % next cell downstream
        flow_direction_parents{nextCell}(flow_direction_parents{nextCell} == curCell) = [];

        flow_direction(curCell) = previousCell; % set the new flow direction;
        % add new flow direction parents
        if isempty(find(flow_direction_parents{previousCell} == curCell, 1))
            flow_direction_parents{previousCell} = [flow_direction_parents{previousCell}, curCell];
        end
        
        previousCell = curCell; % set for the next iteration
        curCell = nextCell;
    end
end