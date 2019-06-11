function[spilloverElevation, cellOverflowInto, pairs] = findSpilloverPairsHist(pairs, dem)
%     pairs = pairs(any(pits(pairs) ~= pitId, 2), :);
    [~, b, c] = unique(pairs(:, 1:2), 'rows'); % get unique values and indexes
    d = accumarray(c, 1); % get counts of each unique value;
    pairs = pairs(b(d<2), :); % for any duplicates (d >= 2), remove both entries 

    [val, ord] = min(max(dem(pairs(:, 1:2)), [], 2));
    spilloverElevation = val;
    cellOverflowInto = pairs(ord, 3);
end