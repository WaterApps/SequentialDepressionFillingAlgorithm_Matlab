function[flow_accumulation] = flowAccumulationAdjacent(flow_direction)
    indexes = reshape(1: numel(flow_direction), size(flow_direction));
    [numrows, numcols] = size(flow_direction);
    
    flow_accumulation = nan(size(flow_direction));
    flow_accumulation(~isnan(flow_direction)) = 0;
    for i = 1 : numel(flow_accumulation)
        if flow_accumulation(i) == 0 % begin recursion if it hasn't been yet visited
            [r, c] = ind2sub(size(flow_direction), i);
            parentSum = recursiveTraverse(r,c);
        end
    end
 
    function tot = recursiveTraverse(r, c)
        flow_accumulation(r,c) =  flow_accumulation(r,c) + 1; % count yourself
        for x = -1 : 1 % loop through neighboring cells
            for y = -1 : 1
                if x == 0 && y ==0
                    continue; % skip center cell of 3x3 neighborhood
                end
                if r+y > numrows || r+y < 1 || c+x > numcols || c+x < 1
                    continue; % skip neighbors outside the dem range
                end
                curCell = indexes(r, c);
                neighbor = indexes(r+y, c+x);
                if flow_direction(neighbor) == curCell
                    if flow_accumulation(neighbor) == 0
                        parentSum = recursiveTraverse(r+y, c+x);
                    else
                        parentSum = flow_accumulation(neighbor);
                    end
                    flow_accumulation(curCell) =  flow_accumulation(curCell) + parentSum;
                end
            end
        end
        tot = flow_accumulation(r,c);
        return
    end
end