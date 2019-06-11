function[spilloverElevation, cellOverflowInto, pairs] = findSpilloverPairsPits(pairs, pits, dem, pitId)
    pairs = pairs(any(pits(pairs) ~= pitId, 2), :);
%     [~, b, c] = unique(lineSegments{second_pit}(:, 1:2), 'rows'); % get unique values and indexes
%     d = accumarray(c, 1); % get counts of each unique value;
%     lineSegments{second_pit} = lineSegments{second_pit}(b(d<2), :); % for any duplicates (d >= 2), remove both entries 

    [val, ord] = min(max(dem(pairs(:, 1:2)), [], 2));
    spilloverElevation = val;
    cellOverflowInto = pairs(ord, pits(pairs(ord, :)) ~= pitId);
end