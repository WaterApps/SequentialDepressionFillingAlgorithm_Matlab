function[flow_accumulation] = flowAccumulation(flow_direction)
    ii = double(flow_direction(~isnan(flow_direction) & flow_direction > 0));
    jj = find(~isnan(flow_direction) & flow_direction > 0);
    flow_direction_parents = cell(size(flow_direction));
    for i = 1 : numel(ii)
        flow_direction_parents{ii(i)} = [flow_direction_parents{ii(i)}, jj(i)];
    end
    flow_accumulation = nan(size(flow_direction));
    flow_accumulation(~isnan(flow_direction)) = 0;
    for i = 1 : length(ii)
        if flow_accumulation(ii(i)) == 0 % begin recursion if it hasn't been yet visited
            parentSum = recursionThree(ii(i));
        end
    end

    function tot = recursionThree(i)
        flow_accumulation(i) =  flow_accumulation(i) + 1; % count yourself
        parents = flow_direction_parents{i};
        for parent = parents
            if (flow_accumulation(parent) == 0)
                parentSum = recursionThree(parent);
            else
                parentSum = flow_accumulation(parent);
            end
            flow_accumulation(i) =  flow_accumulation(i) + parentSum;
        end
        tot = flow_accumulation(i);
        return
    end

%Generate a decent figure
figure(1);
imagesc(flow_accumulation, 'AlphaData',~isnan(flow_accumulation));
title('Flow Accumulation');
cmap = colormap('gray'); % base it on grayscale
colormap(cmap);
h = colorbar;
set(gca,'colorscale','log');



% cutoff = prctile(reshape(temp, [1, numel(temp)]), 99);
% num = 30;
% cmap(1: num+1, :) = repmat([0 : (1/num) : 1.0]', 1, 3);
% cmap(num+2: end) = 1.0;
% ylabel(h, 'foo') % label the colorbar if need be

end
