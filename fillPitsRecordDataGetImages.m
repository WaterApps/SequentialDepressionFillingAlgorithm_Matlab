function[dem, flow_direction, flow_direction_parents, pits, fill_rainfall_amounts, runoff_areas, number_of_pits, mean_depth, std_depth, mean_area, std_area, area] = fillPitsRecordDataGetImages(dem, flow_direction, flow_direction_parents, pits, pitId, pitCell, areaCellCount, spilloverElevation, spilloverTime, volume, filledVolume, outletCell, cellOverflowInto, cellsize, rainfall_intensity, CedarCreekOutletPoint, R, dem_selection)
%First, find the pit containing the Cedar Creek outlet, and set its
%spillover time to infinity.  This will prevent it from overflowing into
%other depressions. Only other depressions may overflow into it.
[cedarCreek, ~] = find(pitId == pits(CedarCreekOutletPoint(1), CedarCreekOutletPoint(2)));
spilloverTime(cedarCreek) = Inf;

% Order pits according to spillover time and initialize the current maximum
% and minimum pit ID numbers
[~, order] = sort(spilloverTime);
pitId = pitId(order);
pitCell = pitCell(order, :);
areaCellCount = areaCellCount(order);
spilloverElevation = spilloverElevation(order);
spilloverTime = spilloverTime(order);
volume = volume(order);
filledVolume = filledVolume(order);
outletCell = outletCell(order, :);
cellOverflowInto = cellOverflowInto(order, :);
max_id = max(pitId);
nanCount = nansum(nansum(~isnan(dem)));

% image_frames = [1, 14275, 14533, 14803, 15034, 15150, 15225]; % 30m - these are
% the chosen frames based on large fill events
% image_frames = [1, 39027, 41217, 43307, 44922, 45653, 45912]; % 10m - these are
% the chosen frames based on large fill events

% image_frames = [1, 14722, 15048, 15149, 15215, 15236]; % 30m - these are
% chosen based on sensible rainfall amounts - 0, 25, 75, 150, 300, 500 (mm)
% image_frames = [1, 42722, 45024, 45647, 45898, 45939]; % 10m - these are
% chosen based on sensible rainfall amounts - 0, 25, 75, 150, 300, 500 (mm)
frames = [1, 14722, 15048, 15149, 15215, 15236;
          1, 42722, 45024, 45647, 45898, 45939;
          1, 443920, 445927, 446509, 446779, 446842];
image_frames = frames(dem_selection, :);
dem_options = ['30m_'; '10m_'; '03m_'];

% Find the total potential mergers and preallocate matrices accordingly.
potential_merges = max_id;
fill_rainfall_amounts = zeros(potential_merges, 1);
runoff_areas = zeros(potential_merges, 1);
number_of_pits = zeros(potential_merges, 1);
mean_depth = zeros(potential_merges, 1);
std_depth = zeros(potential_merges, 1);
mean_area = zeros(potential_merges, 1);
std_area = zeros(potential_merges, 1);
area = zeros(potential_merges, 1);

% Get result values before rainfall begins
fill_rainfall_amounts(1) = 0;
number_of_pits(1) = sum(~isnan(spilloverTime));
runoff_areas(1) = areaCellCount(cedarCreek).*100./nanCount;
mean_depth(1) = nanmean(spilloverTime)*rainfall_intensity;
std_depth(1) = nanstd(spilloverTime)*rainfall_intensity;
mean_area(1) = nanmean(areaCellCount)*0.000247105*cellsize^2; % in acres
std_area(1) = nanstd(areaCellCount)*0.000247105*cellsize^2; % in acres
area(1) = areaCellCount(1)*0.000247105*cellsize^2; % in acres

% Make a movie of the pits and a gif
% puddle_dem = zeros(size(dem,1), size(dem,2), 3);
% filename = '30m_sample_pits.gif';
% dem_filename = '30m_sample_dem.gif';

% pitsMovie(potential_merges) = struct('cdata',[],'colormap',[]);
% % This colormap accounts for negative pitst
% pitsColormap = zeros(abs(min(min(pits))) + max(max(pits)) - 1, 3);
% pitsColormap(abs(min(min(pits))): end, :) = rand(max(max(pits)), 3);
% imagesc(pits)
% axis equal;
% colormap(pitsColormap);
% set(gca,'visible','off')
% set(gca,'position',[0 0 1 1], 'units', 'normalized')
% drawnow
% pitsMovie(1) = getframe;
% im = frame2im(pitsMovie(1));
% [imind,cm] = rgb2ind(im,256);
% imwrite(imind,cm,filename,'gif', 'Loopcount',inf, 'DelayTime', 0.2);
% 
% demMovie(potential_merges-519) = struct('cdata',[],'colormap',[]);
% imagesc(dem)
% axis equal;
% colormap gray;
% set(gca,'visible','off')
% set(gca,'position',[0 0 1 1], 'units', 'normalized')
% drawnow
% demMovie(1) = getframe;
% im = frame2im(demMovie(1));
% [imind,cm] = rgb2ind(im,256);
% imwrite(imind,cm,dem_filename,'gif', 'Loopcount',inf, 'DelayTime', 0.2);
% 
% 
% frame_idx = 2;
pits = uint32(pits);

pitsColormap = rand(max(max(pits))+1, 3); % random color for every depression, plus one for depression ID 0.
pitsColormap(1,:) = 1; % make pit ID 0 white
RGB = ind2rgb(pits, pitsColormap);
image(RGB);
%imagesc(pits)
axis equal;
set(gca,'visible','off');
set(gca,'position',[0 0 1 1], 'units', 'normalized');
drawnow;
imwrite(RGB, strcat(dem_options(dem_selection, :),num2str(1),'.png'));
geotiffwrite(strcat(dem_options(dem_selection, :), num2str(1), '.tif'), pits, R, 'CoordRefSysCode', 26916);

idx = 2;
while ~isnan(spilloverTime(2))
    fill_rainfall_amounts(idx) = spilloverTime(1)*rainfall_intensity;
    area(idx) = areaCellCount(1)*0.000247105*cellsize^2; % in acres
    second_pit_ID = pits(cellOverflowInto(1, 1), cellOverflowInto(1, 2));
    [numrows,numcols] = size(dem);
    [second_pit, ~] = find(pitId == second_pit_ID);
    
    spilloverElevation(second_pit) = NaN;
    outletCell(second_pit, :) = [0, 0];
    cellOverflowInto(second_pit, :) = [0, 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isinf(spilloverTime(second_pit))
        % re-ID first pit, raise/fill first pit cells
        raisedPointsCount = 0;
        indicesToCheck = int16(zeros(areaCellCount(1), 2));
        indicesToCheck(1, :) = pitCell(1, :);
        j = 2;
        for i = 1 : size(indicesToCheck, 1)
            r = indicesToCheck(i, 1);
            c = indicesToCheck(i, 2);
    %         puddle_dem(r,c, :) = [0, 0, 1];
            if (dem(r,c) <= spilloverElevation(1))
                raisedPointsCount = raisedPointsCount + 1;
                dem(r,c) = spilloverElevation(1);
            else
                pits(r,c) = pitId(second_pit);
            end

            if (dem(r,c) < spilloverElevation(second_pit))
                volume(second_pit) = volume(second_pit) + ((spilloverElevation(second_pit) - dem(r,c))*cellsize*cellsize);
            end

            for x = 1 : 3 
                for y = 1 : 3
                    if flow_direction_parents(r,c,y,x)
                        indicesToCheck(j, :) = [r+y-2, c+x-2];
                        j = j + 1;
                    end
                end
            end
        end
        
        % Resolve the flow direction of the filled cells.  Begin by directing
        % the spillover point to the pre-identified outletSpillover point
        indicesToCheck = int16(zeros(raisedPointsCount, 2));
        flow_direction(outletCell(1,1), outletCell(1,2), :) = cellOverflowInto(1,:);
        rr = outletCell(1, 1);
        cc = outletCell(1, 2);
        r = cellOverflowInto(1, 1);
        c = cellOverflowInto(1, 2);
        flow_direction_parents(r,c, rr-r+2, cc-c+2) = true;
        pits(outletCell(1, 1), outletCell(1, 2)) = pitId(second_pit);
        indicesToCheck(1, :) = outletCell(1,:);
        j = 2;
        for i = 1 : size(indicesToCheck, 1)
            r = indicesToCheck(i, 1);
            c = indicesToCheck(i, 2);
            for x = -1 : 1 % loop through neighboring cells
                for y = -1 : 1
                    if x == 0 && y ==0 
                        continue; % skip center cell of 3x3 neighborhood
                    end
                    % Check if the point needs to be resolved, but hasn't
                    % already be resolved. This avoids redundant additions of
                    % points to the "next up" list of points to be resolved.
                    if ((pits(r+y, c+x) == pitId(1)) && (pits(r+y, c+x) ~= pitId(second_pit)))
    %                     if flow_direction(r+y, c+x, :) ~= -1
    %                         flow_direction_parents(flow_direction(r+y, c+x, 1), flow_direction(r+y, c+x, 2), (r+y)-(flow_direction(r+y, c+x, 1))+2, (c+x)-(flow_direction(r+y, c+x, 2))+2) = false;
    %                     end
                        flow_direction(r+y, c+x, :) = [r,c];
    %                     flow_direction_parents(r,c,y+2, x+2) = true;
                        pits(r+y, c+x) = pitId(second_pit);
                        indicesToCheck(j, :) = [r+y, c+x];
                        j = j + 1;
                    end
                end
            end
        end

        % Correct flow direction parents
        for i = 1 : size(indicesToCheck, 1)
            r = indicesToCheck(i, 1);
            c = indicesToCheck(i, 2);
            flow_direction_parents(r,c,:,:) = false;
            for x = -1 : 1 % loop through neighboring cells
                for y = -1 : 1
                    if x == 0 && y == 0 % skip center (target) cell of 3x3 neighborhood
                        continue;
                    end

                    if r+y > numrows || r+y < 1 || c+x > numcols || c+x < 1
                        continue; % skip neighbors outside the matrix range
                    end
                    if flow_direction(r+y, c+x, :) ~= 0
                        if isequal(flow_direction(r+y, c+x, 1), r) && isequal(flow_direction(r+y, c+x, 2), c)
                            flow_direction_parents(r, c, y+2, x+2) = true;
                        end
                    end
                end
            end
        end

        areaCellCount(second_pit) = areaCellCount(second_pit) +  areaCellCount(1);
        filledVolume(second_pit) = volume(1) + filledVolume(second_pit);
        
    else
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        % Begin identifying the spillover location in the first pit.
        indicesToCheck = int16(zeros(areaCellCount(1), 2));
        indicesToCheck(1, :) = pitCell(1, :);
        j = 2;
        for i = 1 : size(indicesToCheck, 1)
            r = indicesToCheck(i, 1);
            c = indicesToCheck(i, 2);
            for x = -1 : 1 % loop through neighboring cells
                for y = -1 : 1
                    if x == 0 && y ==0 
                        continue; % skip center cell of 3x3 neighborhood
                    end
                    if r+y > numrows || r+y < 1 || c+x > numcols || c+x < 1
                        continue; % skip neighbors outside the dem range
                    end
                    if isnan(dem(r+y, c+x))
                        continue;
                    end

                    % First, check for neighbors to check next.
                    if flow_direction_parents(r,c,y+2, x+2)
                        indicesToCheck(j, :) = [r+y, c+x];
                        j = j + 1;
                    end

                    % Next, check if it is on the border.
                    if ((pits(r+y, c+x) ~= pitId(second_pit)) && (pits(r+y, c+x) ~= pitId(1)))
                        % if minimum ridge elevation value is still NaN from
                        % initialization or if the ridge elevations (in and the
                        % neighbor just out of the pit) are BOTH less than the current
                        % minimum ridge elevation in the pit
                        cur_cell_elev = dem(r, c);
                        neighbor_elev = dem(r+y, c+x);
                        if isnan(spilloverElevation(second_pit)) || (cur_cell_elev <= spilloverElevation(second_pit) && neighbor_elev <= spilloverElevation(second_pit))
                            spilloverElevation(second_pit) = max([neighbor_elev, cur_cell_elev]);
                            outletCell(second_pit, :) = [r, c];
                            cellOverflowInto(second_pit, :) = [r+y, c+x];
                        end
                    end
                end
            end
        end

        % Extend the search for the spillover locaiton to the second pit. Only 
        % at the end of this step is the spillover location is known.
        indicesToCheck = int16(zeros(areaCellCount(second_pit), 2));
        indicesToCheck(1, :) = pitCell(second_pit, :);
        j = 2;
        for i = 1 : size(indicesToCheck, 1)
            r = indicesToCheck(i, 1);
            c = indicesToCheck(i, 2);
            for x = -1 : 1 % loop through neighboring cells
                for y = -1 : 1
                    if x == 0 && y ==0 
                        continue; % skip center cell of 3x3 neighborhood
                    end
                    if r+y > numrows || r+y < 1 || c+x > numcols || c+x < 1
                        continue; % skip neighbors outside the dem range
                    end
                    if isnan(dem(r+y, c+x))
                        continue;
                    end

                    % First, check for neighbors to check next.
                    if flow_direction_parents(r,c,y+2, x+2)
                        indicesToCheck(j, :) = [r+y, c+x];
                        j = j + 1;
                    end

                    % Next, check if it is on the border.
                    if ((pits(r+y, c+x) ~= pitId(second_pit)) && (pits(r+y, c+x) ~= pitId(1)))
                        % if minimum ridge elevation value is still NaN from
                        % initialization or if the ridge elevations (in and the
                        % neighbor just out of the pit) are BOTH less than the current
                        % minimum ridge elevation in the pit
                        cur_cell_elev = dem(r, c);
                        neighbor_elev = dem(r+y, c+x);
                        if isnan(spilloverElevation(second_pit)) || (cur_cell_elev <= spilloverElevation(second_pit) && neighbor_elev <= spilloverElevation(second_pit))
                            spilloverElevation(second_pit) = max([neighbor_elev, cur_cell_elev]);
                            outletCell(second_pit, :) = [r, c];
                            cellOverflowInto(second_pit, :) = [r+y, c+x];
                        end
                    end
                end
            end
        end

        % Update the filled volumes of the pit
        filledVolume(second_pit) = volume(1) + filledVolume(second_pit);
        volume(second_pit) = filledVolume(second_pit);

        % re-ID first pit, raise/fill first pit cells
        raisedPointsCount = 0;
        indicesToCheck = int16(zeros(areaCellCount(1), 2));
        indicesToCheck(1, :) = pitCell(1, :);
        j = 2;
        for i = 1 : size(indicesToCheck, 1)
            r = indicesToCheck(i, 1);
            c = indicesToCheck(i, 2);
    %         puddle_dem(r,c, :) = [0, 0, 1];
            if (dem(r,c) <= spilloverElevation(1))
                raisedPointsCount = raisedPointsCount + 1;
                dem(r,c) = spilloverElevation(1);
            else
                pits(r,c) = pitId(second_pit);
            end

            if (dem(r,c) < spilloverElevation(second_pit))
                volume(second_pit) = volume(second_pit) + ((spilloverElevation(second_pit) - dem(r,c))*cellsize*cellsize);
            end

            for x = 1 : 3 
                for y = 1 : 3
                    if flow_direction_parents(r,c,y,x)
                        indicesToCheck(j, :) = [r+y-2, c+x-2];
                        j = j + 1;
                    end
                end
            end
        end

        % Calculate volume through the second pit
        indicesToCheck = int16(zeros(areaCellCount(second_pit), 2));
        indicesToCheck(1, :) = pitCell(second_pit, :);
        j = 2;
        for i = 1 : size(indicesToCheck, 1)
            r = indicesToCheck(i, 1);
            c = indicesToCheck(i, 2);
            if (dem(r,c) < spilloverElevation(second_pit))
                volume(second_pit) = volume(second_pit) + ((spilloverElevation(second_pit) - dem(r,c))*cellsize*cellsize);
            end

            for x = 1 : 3 
                for y = 1 : 3
                    if flow_direction_parents(r,c,y,x)
                        indicesToCheck(j, :) = [r+y-2, c+x-2];
                        j = j + 1;
                    end
                end
            end
        end

        % Resolve the flow direction of the filled cells.  Begin by directing
        % the spillover point to the pre-identified outletSpillover point
        indicesToCheck = int16(zeros(raisedPointsCount, 2));
    %     if flow_direction(outletCell(1,1), outletCell(1,2), :) ~= -1
    %         flow_direction_parents(flow_direction(outletCell(1,1), outletCell(1,2), 1), flow_direction(outletCell(1,1), outletCell(1,2), 2), (outletCell(1,1))-(flow_direction(outletCell(1,1), outletCell(1,1), 1))+2, (outletCell(1,2))-(flow_direction(outletCell(1,1), outletCell(1,2), 2))+2) = false;
    %     end
        flow_direction(outletCell(1,1), outletCell(1,2), :) = cellOverflowInto(1,:);
        rr = outletCell(1, 1);
        cc = outletCell(1, 2);
        r = cellOverflowInto(1, 1);
        c = cellOverflowInto(1, 2);
        flow_direction_parents(r,c, rr-r+2, cc-c+2) = true;
        pits(outletCell(1, 1), outletCell(1, 2)) = pitId(second_pit);
        indicesToCheck(1, :) = outletCell(1,:);
        j = 2;
        for i = 1 : size(indicesToCheck, 1)
            r = indicesToCheck(i, 1);
            c = indicesToCheck(i, 2);
            for x = -1 : 1 % loop through neighboring cells
                for y = -1 : 1
                    if x == 0 && y ==0 
                        continue; % skip center cell of 3x3 neighborhood
                    end
                    % Check if the point needs to be resolved, but hasn't
                    % already be resolved. This avoids redundant additions of
                    % points to the "next up" list of points to be resolved.
                    if ((pits(r+y, c+x) == pitId(1)) && (pits(r+y, c+x) ~= pitId(second_pit)))
    %                     if flow_direction(r+y, c+x, :) ~= -1
    %                         flow_direction_parents(flow_direction(r+y, c+x, 1), flow_direction(r+y, c+x, 2), (r+y)-(flow_direction(r+y, c+x, 1))+2, (c+x)-(flow_direction(r+y, c+x, 2))+2) = false;
    %                     end
                        flow_direction(r+y, c+x, :) = [r,c];
    %                     flow_direction_parents(r,c,y+2, x+2) = true;
                        pits(r+y, c+x) = pitId(second_pit);
                        indicesToCheck(j, :) = [r+y, c+x];
                        j = j + 1;
                    end
                end
            end
        end

        % Correct flow direction parents
        for i = 1 : size(indicesToCheck, 1)
            r = indicesToCheck(i, 1);
            c = indicesToCheck(i, 2);
            flow_direction_parents(r,c,:,:) = false;
            for x = -1 : 1 % loop through neighboring cells
                for y = -1 : 1
                    if x == 0 && y == 0 % skip center (target) cell of 3x3 neighborhood
                        continue;
                    end

                    if r+y > numrows || r+y < 1 || c+x > numcols || c+x < 1
                        continue; % skip neighbors outside the matrix range
                    end
                    if flow_direction(r+y, c+x, :) ~= 0
                        if isequal(flow_direction(r+y, c+x, 1), r) && isequal(flow_direction(r+y, c+x, 2), c)
                            flow_direction_parents(r, c, y+2, x+2) = true;
                        end
                    end
                end
            end
        end

        areaCellCount(second_pit) = areaCellCount(second_pit) +  areaCellCount(1);
        spilloverTime(second_pit) = volume(second_pit)/(rainfall_intensity*areaCellCount(second_pit)*(cellsize^2));

        if spilloverTime(second_pit) < 0
            spilloverTime(second_pit) = Inf;
        end
    end
    
    spilloverTime(1) = NaN;
    [~, order] = sort(spilloverTime);
    pitId = pitId(order);
    pitCell = pitCell(order, :);
    areaCellCount = areaCellCount(order);
    spilloverElevation = spilloverElevation(order);
    spilloverTime = spilloverTime(order);
    volume = volume(order);
    filledVolume = filledVolume(order);
    outletCell = outletCell(order, :);
    cellOverflowInto = cellOverflowInto(order, :);
    
%     if idx >=519
%         imagesc(pits)
%         axis equal;
%         colormap(pitsColormap);
%         set(gca,'visible','off')
%         set(gca,'position',[0 0 1 1], 'units', 'normalized')
%         drawnow
%         pitsMovie(frame_idx) = getframe;
%         im = frame2im(pitsMovie(frame_idx));
%         [imind,cm] = rgb2ind(im,256);
%         imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime', 0.2);
% 
%         imagesc(dem)
%         axis equal;
%         colormap(gray);
%         set(gca,'visible','off')
%         set(gca,'position',[0 0 1 1], 'units', 'normalized')
%         drawnow
%         demMovie(frame_idx) = getframe;
%         im = frame2im(demMovie(frame_idx));
%         [imind,cm] = rgb2ind(im,256);
%         imwrite(imind,cm,dem_filename,'gif','WriteMode','append', 'DelayTime', 0.2);
%         
%         frame_idx = frame_idx + 1;
%     end
    
    %%
    [cedarCreek, ~] = find(pitId == pits(CedarCreekOutletPoint(1), CedarCreekOutletPoint(2)));
    number_of_pits(idx) = sum(~isnan(spilloverTime));
    runoff_areas(idx) = areaCellCount(cedarCreek).*100./nanCount;
    mean_depth(idx) = nanmean(spilloverTime)*rainfall_intensity;
    std_depth(idx) = nanstd(spilloverTime)*rainfall_intensity;
    mean_area(idx) = nanmean(areaCellCount)*0.000247105*cellsize^2;
    std_area(idx) = nanstd(areaCellCount)*0.000247105*cellsize^2;

    if (any(idx == image_frames))
        RGB = ind2rgb(pits, pitsColormap);
        image(RGB);
        axis equal;
        set(gca,'visible','off');
        set(gca,'position',[0 0 1 1], 'units', 'normalized');
        drawnow;
        imwrite(RGB, strcat(dem_options(dem_selection, :), num2str(idx),'.png'));
        geotiffwrite(strcat(dem_options(dem_selection, :), num2str(idx), '.tif'), pits, R, 'CoordRefSysCode', 26916);
    end
    
    % This stops execution after the last frame is gathered to save time.
    if (idx == image_frames(end))
        disp('Done');
        break;
    end
    idx = idx + 1
end

% movie(pitsMovie, 5);
% movie(demMovie, 5);

% save pitsMovie;
% save demMovie;

end

% rows = (r-1)*(sizGGetImageGe(imagery,1)/size(puddle_dem,1))+1:r*(size(imagery,1)/size(puddle_dem,1));
% cols = (c-1)*(size(imagery,1)/size(puddle_dem,1))+1:c*(size(imagery,1)/size(puddle_dem,1));
% puddle_dem(rows, cols) 




