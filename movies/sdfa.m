function[fill_dem, fill_flow_direction, fill_pits, fill_flow_accumulation, depthFlow, rainfall_excess, runoff] = sdfa(filepath, fillRainfallExcess);

  defaultFilePath = './Feldun.tif';
  defaultFillAmount = Inf;

  if ~exist('filepath', 'var')
    filepath = defaultFilePath;
  end

  if ~exist('fillRainfallExcess', 'var')
    fillRainfallExcess = Inf;
  end

  visualize_merging = true; % rendering slows the filling process

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
  [pits, pairs, cellIndexes, pitId, pitCell, edgePitCell, areaCellCount, spilloverElevation, vca, volume, filledVolume, cellOverflowInto, perimeterCell] = ...
      getDepressions(dem, flow_direction, flow_direction_parents, R.CellExtentInWorldX);

  %% Fill Pits
  [fill_dem, fill_flow_direction, fill_pits, depthFlow, rainfall_excess, runoff] = ...
      fillDepressionsMovieFlowAccumulation(fillRainfallExcess, dem, flow_direction, pits, pairs, cellIndexes, pitId, pitCell, edgePitCell, areaCellCount, spilloverElevation, vca, volume, filledVolume, cellOverflowInto, perimeterCell, R, visualize_merging);

  %% Flow Accumulation
  disp('Computing Flow Accumulation')
  [fill_flow_accumulation] = flowAccumulation(fill_flow_direction);

  [fpath, name, ext] = fileparts(filepath)
  ginfo = geotiffinfo(filepath);
  tag = ginfo.GeoTIFFTags.GeoKeyDirectoryTag
  if (isfield(tag, 'Unknown'))
    tag = rmfield(tag, 'Unknown')
  end
  geotiffwrite(strcat(fpath, filesep, name, '_flowAccumulation.tif'), fill_flow_accumulation, R, 'GeoKeyDirectoryTag', tag);
  geotiffwrite(strcat(fpath, filesep, name, '_fill_dem.tif'), fill_dem, R, 'GeoKeyDirectoryTag', tag);
  geotiffwrite(strcat(fpath, filesep, name, '_catchments.tif'), fill_pits, R, 'GeoKeyDirectoryTag', tag);

end