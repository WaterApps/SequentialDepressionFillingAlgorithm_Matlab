function[fill_dem, fill_flow_direction, fill_pits, fill_flow_accumulation] = sdfa(filepath, fillRainfallExcess);

  defaultFilePath = './Feldun.tif';
  defaultFillAmount = Inf;

  if ~exist('filepath', 'var')
    filepath = defaultFilePath;
  end

  if ~exist('fillRainfallExcess', 'var')
    fillRainfallExcess = Inf;
  end

  visualize_merging = false; % rendering slows the filling process

  %%
  [dem, R] = geotiffread(filepath);
  cellsize = R.CellExtentInWorldX; % DEM cellsize in meters
  disp(['DEM loaded...dimensions are: ', num2str(size(dem))]);

  %% The clipped DEMs have an assigned NaN value of -3.4028e38. Find and mark as matlab NaN.
  dem = double(dem);
  parfor i = 1 : numel(dem)
      if (dem(i) < -100)
          dem(i) = NaN
      end
  end

  %% Recompute Flow Direction
  disp('Flow direction matrix was NOT found')
  disp('Computing Flow Direction');
  [flow_direction] = d8FlowDirection(dem);
  [flow_direction_parents] = d8FlowDirectionParents(flow_direction);

  %% Identify Pits, Compute Matrix/Map with Pit ID for each cell
  disp('Creating Pit Dataset');
  [pits, pairs, cellIndexes, pitId, pitCell, areaCellCount, spilloverElevation, vca, volume, filledVolume, cellOverflowInto] = ...
      getDepressions(dem, flow_direction, flow_direction_parents, R.CellExtentInWorldX);

  %% Fill Pits
  [fill_dem, fill_flow_direction, fill_pits] = ...
      fillDepressions(fillRainfallExcess, dem, flow_direction, pits, pairs, cellIndexes, pitId, pitCell, areaCellCount, spilloverElevation, vca, volume, filledVolume, cellOverflowInto, R, visualize_merging);

  %% Flow Accumulation
  disp('Computing Flow Accumulation')
  [fill_flow_accumulation] = flowAccumulation(fill_flow_direction);

  [filepath, name, ext] = fileparts(filepath)
  geotiffwrite(strcat(filepath,name, '_flowAccumulation.tif'), flow_accumulation, R);
  geotiffwrite(strcat(filepath,name, '_fill_dem.tif'), fill_dem, R);
  geotiffwrite(strcat(filepath,name, '_catchments.tif'), fill_pits, R);

end
