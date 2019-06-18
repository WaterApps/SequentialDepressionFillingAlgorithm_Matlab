function[flow_direction] = updateFlowDirectionAdjacentDirect(cellOverflowInto, perimeterCell, flow_direction)
    % reorient filled zone back toward spillover cell, change pit ID to
    % second_pit, and raise DEM
    previousCell = cellOverflowInto;
    curCell = perimeterCell;
    while flow_direction(curCell) > 0
        nextCell = flow_direction(curCell); % store the old flow direction;
        flow_direction(curCell) = previousCell; % set the new flow direction;
        previousCell = curCell; % advance for the next iteration
        curCell = nextCell; % advance for the next iteration
    end
end