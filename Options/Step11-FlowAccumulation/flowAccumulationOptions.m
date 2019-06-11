function[flow_accumulation] = flowAccumulationOptions(flow_direction)

% [fa1] = flowAccumulationSparseParents(flow_direction);
[fa2] = flowAccumulationDownstreamIncrement(flow_direction);
[fa3] = flowAccumulationDownstreamIncrement2(flow_direction);
[flow_accumulation] = flowAccumulationCellArrayParents(flow_direction);
[fa4] = flowAccumulationAdjacent(flow_direction);

if ~isequaln(fa2, fa3, flow_accumulation, fa4)
    
end

%Generate a decent figure
% temp = flow_accumulation;
% figure(1);
% imagesc(temp, 'AlphaData',~isnan(temp));
% title('Flow Accumulation');
% cmap = colormap('gray'); % base it on grayscale
% colormap(cmap);
% h = colorbar;
% set(gca,'colorscale','log');
% cutoff = prctile(reshape(temp, [1, numel(temp)]), 99);
% num = 30;
% cmap(1: num+1, :) = repmat([0 : (1/num) : 1.0]', 1, 3);
% cmap(num+2: end) = 1.0;
% ylabel(h, 'foo') % label the colorbar if need be

end
