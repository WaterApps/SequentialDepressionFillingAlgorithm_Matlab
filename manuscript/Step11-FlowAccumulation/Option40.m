function[flow_accumulation] = flowAccumulationSparseParents(flow_direction)
ii = double(flow_direction(~isnan(flow_direction) & flow_direction > 0)); % the cells that flow elsewhere
jj = find(~isnan(flow_direction) & flow_direction > 0); % the destinations of the cells from ii
flow_direction_parents = sparse(ii, jj, true, numel(flow_direction), numel(flow_direction), length(jj));
flow_accumulation = nan(size(flow_direction));
flow_accumulation(~isnan(flow_direction)) = 0;
for i = 1 : length(ii)
    if flow_accumulation(ii(i)) == 0 % begin recursion if it hasn't been yet visited
        parentSum = recursionTwo(ii(i));
    end
end

function tot = recursionTwo(curCell)
    flow_accumulation(curCell) =  flow_accumulation(curCell) + 1; % count yourself6.8m/81.2m/7.3m
    parents = [find(flow_direction_parents(curCell, :))];
    for parent = parents
        if (flow_accumulation(parent) == 0)
            parentSum = recursionTwo(parent);
        else
            parentSum = flow_accumulation(parent);
        end
        flow_accumulation(curCell) =  flow_accumulation(curCell) + parentSum;
    end
    tot = flow_accumulation(curCell);
    return %tot;
end
end