function [flow_direction_parents] = d8FlowDirectionParentsCellArray(flow_direction)

%Now find parents
fprintf('Creating flow direction parents adjacency matrix...')
flow_direction_parents = cell(size(flow_direction));
for i = 1 : numel(flow_direction)
    if isnan(flow_direction(i))
       continue; 
    end
    if flow_direction(i) > 0
        child = flow_direction(i);
        flow_direction_parents{child} = [flow_direction_parents{child}, i];
    end
end
fprintf('Done')
fprintf('\n')

end

