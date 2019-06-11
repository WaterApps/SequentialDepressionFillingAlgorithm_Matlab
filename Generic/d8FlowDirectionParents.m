function [flow_direction_parents] = d8FlowDirectionParents(flow_direction)

%Now find parents
fprintf('Creating flow direction parents adjacency matrix...')
ii = double(flow_direction(~isnan(flow_direction) & flow_direction > 0)); % the cells that flow elsewhere
jj = find(~isnan(flow_direction) & flow_direction > 0); % the destinations of the cells from ii
flow_direction_parents = cell(size(flow_direction));
for i = 1 : numel(ii)
    flow_direction_parents{ii(i)} = [flow_direction_parents{ii(i)}, jj(i)];
end
fprintf('Done')
fprintf('\n')

end

