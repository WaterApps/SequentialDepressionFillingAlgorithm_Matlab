function[flow_direction] = updateFlowDirectionNonAdjacent(cellOverflowInto, pitCell, flow_direction)
    flow_direction(pitCell) = cellOverflowInto;
end