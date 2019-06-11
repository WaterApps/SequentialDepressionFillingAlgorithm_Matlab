function[flow_accumulation] = flowAccumulationDownstreamIncrement2(flow_direction)
    ii = double(flow_direction(~isnan(flow_direction) & flow_direction > 0));
    jj = find(~isnan(flow_direction) & flow_direction > 0);
    flow_accumulation = nan(size(flow_direction));
    flow_accumulation(~isnan(flow_direction)) = 1;
    
    flow_direction_parents = cell(size(flow_direction));
    for i = 1 : numel(ii)
        flow_direction_parents{ii(i)} = [flow_direction_parents{ii(i)}, jj(i)];
    end

    ii = find(~isnan(flow_direction) & ~cellfun(@isempty,flow_direction_parents));
    for i = 1 : length(ii)
    	curCell = ii(i);
        while flow_direction(curCell) > 0
            curCell = flow_direction(curCell);
            flow_accumulation(curCell) = flow_accumulation(curCell) +1;
        end
    end
end