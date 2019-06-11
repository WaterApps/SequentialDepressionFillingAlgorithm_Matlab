function[flow_accumulation] = flowAccumulationDownstreamIncrement(flow_direction)
    flow_accumulation = zeros(size(flow_direction));

    for i = 1 : numel(flow_direction)
        if ~isnan(flow_direction(i))
            flow_accumulation(i) = flow_accumulation(i) + 1;
            curCell = i;
            while flow_direction(curCell) > 0
                curCell = flow_direction(curCell);
                flow_accumulation(curCell) = flow_accumulation(curCell) +1; % increment one for the parent flowing to it (why we recursed here)
            end
        end
    end


end