% thirty = load('CedarUpper_30m');
% ten = load('CedarUpper_10m');
% three = load('CedarUpper_03m');
    
close all;
sum(~isnan(thirty.dem(:)))./length(thirty.rainfall_excess)
sum(~isnan(ten.dem(:)))./length(ten.rainfall_excess)
sum(~isnan(three.dem(:)))./length(three.rainfall_excess)

% Number of Depressions (absolute count) 
figure(1);
scatter(thirty.rainfall_excess.*1000, thirty.number_of_pits./1000, '.', 'red');
hold on;
scatter(ten.rainfall_excess.*1000, ten.number_of_pits./1000, '.', 'blue');
scatter(three.rainfall_excess.*1000, three.number_of_pits./1000, '.', 'green');
title('Number of Depressions Remaining');
legend('30m DEM', '10m DEM', '3m DEM');
xlabel('Rainfall Excess (mm)');
ylabel('Number of Depressions Remaining (1000s)');
xlim([0,100]);
ylim([0,700]);

% Percent Running Off
% Convert to percentage
thirty.percent_running_off = thirty.runoff_areas.*100./thirty.runoff_areas(end);
ten.percent_running_off = ten.runoff_areas.*100./ten.runoff_areas(end);
three.percent_running_off = three.runoff_areas.*100./three.runoff_areas(end);
figure(2);
stairs(thirty.rainfall_excess.*1000, thirty.percent_running_off, '-r', 'LineWidth', 1.5);
hold on;
stairs(ten.rainfall_excess.*1000, ten.percent_running_off, '--b', 'LineWidth', 1.5);
stairs(three.rainfall_excess.*1000, three.percent_running_off , ':g', 'LineWidth', 1.5);
title('Percent of Watershed Draining to Outlet');
legend('30m DEM', '10m DEM', '3m DEM', 'Location', 'southeast');
ylabel('Percent of Watershed Area Running Off (%)');
xlabel('Rainfall Excess (mm)');
%Going to need some limits to the plot
xlim([0, 500]);
ylim([0, 100]);
% a = [0, 25, 75, 150, 300, 500];
% for i = 1 : length(a)
%     plot([a(i), a(i)], [0, 110]);
% end

% Storage (m3)
figure(3);
stairs(thirty.rainfall_excess.*1000, thirty.storage_volume, '-r', 'LineWidth', 1.5);
hold on;
stairs(ten.rainfall_excess.*1000, ten.storage_volume, '--b', 'LineWidth', 1.5);
stairs(three.rainfall_excess.*1000, three.storage_volume , ':g', 'LineWidth', 1.5);
title('Depression Storage');
ylabel('Depression Storage (m^3)');
xlabel('Rainfall Excess (mm)');
legend('30m DEM', '10m DEM', '3m DEM', 'Location', 'southeast');
xlim([0,500]);
% ylim([0,100000]);

thirty.runoff_volume = thirty.rainfall_excess.*(sum(~isnan(thirty.dem(:)))*30*30) - thirty.storage_volume;
ten.runoff_volume = ten.rainfall_excess.*(sum(~isnan(ten.dem(:)))*10*10) - ten.storage_volume;
three.runoff_volume = three.rainfall_excess.*(sum(~isnan(three.dem(:)))*3*3) - three.storage_volume;
% Storage and Runoff (m3)
figure(4);
stairs(thirty.rainfall_excess.*1000, thirty.storage_volume, '-r', 'LineWidth', 1.5);
hold on;
stairs(ten.rainfall_excess.*1000, ten.storage_volume, '--b', 'LineWidth', 1.5);
stairs(three.rainfall_excess.*1000, three.storage_volume , ':g', 'LineWidth', 1.5);
stairs(thirty.rainfall_excess.*1000, thirty.runoff_volume, '-r', 'LineWidth', 1.5);
stairs(ten.rainfall_excess.*1000, ten.runoff_volume, '--b', 'LineWidth', 1.5);
stairs(three.rainfall_excess.*1000, three.runoff_volume , ':g', 'LineWidth', 1.5);
title('Depression Storage');
ylabel('Depression Storage (m^3)');
xlabel('Rainfall Excess (mm)');
legend('30m DEM', '10m DEM', '3m DEM', 'Location', 'southeast');
xlim([0,500]);
ylim([0,100000000]);

thirty.runoff_volume = thirty.rainfall_excess.*(sum(~isnan(thirty.dem(:)))*30*30) - thirty.storage_volume;
ten.runoff_volume = ten.rainfall_excess.*(sum(~isnan(ten.dem(:)))*10*10) - ten.storage_volume;
three.runoff_volume = three.rainfall_excess.*(sum(~isnan(three.dem(:)))*3*3) - three.storage_volume;
% Storage and Runoff (m3)
figure(5);
stairs(thirty.rainfall_excess.*1000, 1000.*thirty.storage_volume./(sum(~isnan(thirty.dem(:)))*30*30), '-r', 'LineWidth', 1.5);
hold on;
stairs(three.rainfall_excess.*1000, 1000.*three.storage_volume./(sum(~isnan(three.dem(:)))*3*3) , ':g', 'LineWidth', 1.5);
stairs(ten.rainfall_excess.*1000, 1000.*ten.storage_volume./(sum(~isnan(ten.dem(:)))*10*10), '--b', 'LineWidth', 1.5);
% xlim([0,500]);
% ax1 = gca;
% ax1_pos = ax1.Position; % position of first axes
% ax2 = axes('Position', ax1_pos, 'XAxisLocation', 'top', 'YAxisLocation', 'right', 'Color', 'none');
% ax2.XLim = [0 500];
% ax2.YLim = [0 500];
stairs(thirty.rainfall_excess.*1000, 1000.*thirty.runoff_volume./(sum(~isnan(thirty.dem(:)))*30*30), '-r', 'LineWidth', 1.5);
stairs(three.rainfall_excess.*1000, 1000.*three.runoff_volume./(sum(~isnan(three.dem(:)))*3*3), ':g', 'LineWidth', 1.5);
stairs(ten.rainfall_excess.*1000, 1000.*ten.runoff_volume./(sum(~isnan(ten.dem(:)))*10*10), '--b', 'LineWidth', 1.5);
title('Storage Depth and Runoff Depth');
ylabel('Storage Depth and Runoff Depth (mm)');
xlabel('Rainfall Excess (mm)');
legend('30m DEM', '10m DEM', '3m DEM', 'Location', 'northwest');
xlim([0,500]);
% ylim([0,100000000]);

% Plot of storage/runoff (%)
thirty.storage_runoff = thirty.storage_volume.*100./(thirty.storage_volume + thirty.runoff_volume);
ten.storage_runoff = ten.storage_volume.*100./(ten.storage_volume + ten.runoff_volume);
three.storage_runoff = three.storage_volume.*100./(three.storage_volume + three.runoff_volume);

figure(6);
area(ten.rainfall_excess.*1000, ten.storage_runoff, 'FaceColor', 'green');
hold on;
area(three.rainfall_excess.*1000, three.storage_runoff, 'FaceColor', 'blue');
area(thirty.rainfall_excess.*1000, thirty.storage_runoff,  'FaceColor', 'red');
title('Depression Storage vs Runoff');
ylabel('Depression Storage (% of Applied Rainfall Excess)');
xlabel('Rainfall Excess (mm)');
legend('10m DEM', '3m DEM', '30m DEM', 'Location', 'northeast');
xlim([0, 500]);
dim = [.17 .12 .15 .15];
str = 'Depression Storage';
annotation('textbox',dim,'String',str,'FitBoxToText','on', 'EdgeColor', 'none', 'BackgroundColor', 'none', 'Color', 'white', 'FontWeight', 'bold');
dim = [.60 .55 .15 .15];
str = 'Runoff';
annotation('textbox',dim,'String',str,'FitBoxToText','on', 'EdgeColor', 'none', 'BackgroundColor', 'none', 'FontWeight', 'bold');
ylim([0, 100]);

max_value = 300;
figure(7);
area(ten.rainfall_excess.*1000, 1000.*ten.storage_volume./(sum(~isnan(ten.dem(:)))*10*10), 'FaceColor', 'green');
hold on;
area(three.rainfall_excess.*1000, 1000.*three.storage_volume./(sum(~isnan(three.dem(:)))*3*3), 'FaceColor', 'blue');
area(thirty.rainfall_excess.*1000, 1000.*thirty.storage_volume./(sum(~isnan(thirty.dem(:)))*30*30),  'FaceColor', 'red');
title('Depression Storage and Runoff');
ylabel('Depression Storage Depth (mm)');
xlabel('Rainfall Excess (mm)');
legend('10m DEM', ' 3m DEM', '30m DEM', 'Location', 'northwest');
plot([0, max_value], [0, max_value]);
xlim([0, max_value]);
dim = [.5 .17 .15 .15];
str = 'Depression Storage';
annotation('textbox',dim,'String',str,'FitBoxToText','on', 'EdgeColor', 'none', 'BackgroundColor', 'none', 'Color', 'white');
dim = [.7 .47 .15 .15];
str = 'Runoff';
annotation('textbox',dim,'String',str,'FitBoxToText','on', 'EdgeColor', 'none', 'BackgroundColor', 'none');
ylim([0, max_value]);
xlimi=get(gca,'XLim');
ylimi=get(gca,'YLim');
ht = text(0.35*xlimi(2),0.39*ylimi(2),'Total Rainfall Excess');
set(ht,'Rotation',38.5);
set(ht,'FontSize',11);

figure(8);
[nr,nc] = size(thirty.dem);
imagesc(thirty.fill_dem-thirty.dem);
pcolor([thirty.fill_dem-thirty.dem nan(nr,1); nan(1,nc+1)]);
shading flat;
set(gca,'YDir','reverse')
axis equal;
set(gca,'visible','off');
set(gca,'position',[0 0 1 1], 'units', 'normalized');

figure(9);
[nr,nc] = size(ten.dem);
imagesc(ten.fill_dem-ten.dem);
pcolor([ten.fill_dem-ten.dem nan(nr,1); nan(1,nc+1)]);
shading flat;
set(gca,'YDir','reverse')
axis equal;
set(gca,'visible','off');
set(gca,'position',[0 0 1 1], 'units', 'normalized');

figure(10);
[nr,nc] = size(three.dem);
imagesc(three.fill_dem-three.dem);
pcolor([three.fill_dem-three.dem nan(nr,1); nan(1,nc+1)]);
shading flat;
set(gca,'YDir','reverse')
axis equal;
set(gca,'visible','off');
set(gca,'position',[0 0 1 1], 'units', 'normalized');
