function [flow_direction_parents] = d8FlowDirectionParentsSparse(flow_direction)

%Now find parents
fprintf('Creating flow direction parents adjacency matrix...')
ii = double(flow_direction(~isnan(flow_direction) & flow_direction > 0)); % the cells that flow elsewhere
jj = find(~isnan(flow_direction) & flow_direction > 0); % the destinations of the cells from ii
flow_direction_parents = sparse(ii, jj, true, numel(flow_direction), numel(flow_direction), length(ii));
fprintf('Done')
fprintf('\n')

end

