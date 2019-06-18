function[contributingArea] = flowAccumulation(flow_direction_parents, depthFlow)
  contributingArea = false(size(dem));
  coords = = round(ginput(1));
  j = 1; % current length of contributingCells
  i = 1; % iterator
  chunk = 50;
  contributingCells = nan(chunk, 1);
  contributingCells(1) = coords;
  contributingArea(coords) = true;
  while i <= j
      parents = flow_direction_parents{contributingCells(i)};
      contributingArea(parents) = true;
      k = j + length(parents);
      if (k > chunk)
          contributingCells(chunk+1:chunk+50) = zeros(50, 1);
          chunk = chunk + 50;
      end
      contributingCells(j+1 : k) = parents;
      j = k;
      i = i+1;
  end
  contributingCells(j+1:end) = [];


  imagesc(contributingArea);

end
