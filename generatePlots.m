% This file is meant to be ran with outputs already generated for all three
% DEMs. It has been made not to crash if one or more output files are 
% missing, but some legend entries and scales of the plots were intended to
% fit the contents of all three DEMs.

if exist('CedarUpper_30m_outputs.mat', 'file')
    thirty = load('CedarUpper_30m_outputs');
end

if exist('CedarUpper_10m_outputs', 'file')
    ten = load('CedarUpper_10m_outputs');
end

if exist('CedarUpper_03m_outputs', 'file')
    three = load('CedarUpper_03m_outputs');
end

close all;
if exist('thirty', 'var')
    sum(~isnan(thirty.dem(:)))./length(thirty.rainfall_excess)
end
if exist('ten', 'var')
    sum(~isnan(ten.dem(:)))./length(ten.rainfall_excess)
end
if exist('three', 'var')
    sum(~isnan(three.dem(:)))./length(three.rainfall_excess)
end

% Number of Depressions (absolute count) 
figure(1);
if exist('thirty', 'var')
    scatter(thirty.rainfall_excess.*1000, thirty.number_of_pits./1000, '.', 'red');
end
hold on;
if exist('ten', 'var')
    scatter(ten.rainfall_excess.*1000, ten.number_of_pits./1000, '.', 'blue');
end
if exist('three', 'var')
    scatter(three.rainfall_excess.*1000, three.number_of_pits./1000, '.', 'green');
end
title('Number of Depressions Remaining');
legend('30m DEM', '10m DEM', '3m DEM');
xlabel('Rainfall Excess (mm)');
ylabel('Number of Depressions Remaining (1000s)');
xlim([0,100]);
ylim([0,700]);

% Percent Running Off
% Convert to percentage
if exist('thirty', 'var')
    thirty.percent_running_off = thirty.runoff_areas.*100./thirty.runoff_areas(end);
end
if exist('ten', 'var')
    ten.percent_running_off = ten.runoff_areas.*100./ten.runoff_areas(end);
end
if exist('three', 'var')
    three.percent_running_off = three.runoff_areas.*100./three.runoff_areas(end);
end
figure(2);
if exist('thirty', 'var')
    stairs(thirty.rainfall_excess.*1000, thirty.percent_running_off, '-r', 'LineWidth', 1.5);
end
hold on;
if exist('ten', 'var')
    stairs(ten.rainfall_excess.*1000, ten.percent_running_off, '--b', 'LineWidth', 1.5);
end
if exist('three', 'var')
    stairs(three.rainfall_excess.*1000, three.percent_running_off , ':g', 'LineWidth', 1.5);
end
title('Percent of Watershed Draining to Outlet');
legend('30m DEM', '10m DEM', '3m DEM', 'Location', 'southeast');
ylabel('Percent of Watershed Area Running Off (%)');
xlabel('Rainfall Excess (mm)');
%Going to need some limits to the plot
xlim([0, 500]);
ylim([0, 100]);

% Storage (m3)
figure(3);
if exist('thirty', 'var')
    stairs(thirty.rainfall_excess.*1000, thirty.storage_volume, '-r', 'LineWidth', 1.5);
end
hold on;
if exist('ten', 'var')
    stairs(ten.rainfall_excess.*1000, ten.storage_volume, '--b', 'LineWidth', 1.5);
end
if exist('three', 'var')
    stairs(three.rainfall_excess.*1000, three.storage_volume , ':g', 'LineWidth', 1.5);
end
title('Depression Storage');
ylabel('Depression Storage (m^3)');
xlabel('Rainfall Excess (mm)');
legend('30m DEM', '10m DEM', '3m DEM', 'Location', 'southeast');
xlim([0,500]);
% ylim([0,100000]);

% Plot of storage/runoff (%)
if exist('thirty', 'var')
    thirty.storage_runoff = thirty.storage_volume.*100./(thirty.storage_volume + thirty.runoff_volume);
end
if exist('ten', 'var')
    ten.storage_runoff = ten.storage_volume.*100./(ten.storage_volume + ten.runoff_volume);
end
if exist('three', 'var')
    three.storage_runoff = three.storage_volume.*100./(three.storage_volume + three.runoff_volume);
end

figure(6);
if exist('ten', 'var')
    area(ten.rainfall_excess.*1000, ten.storage_runoff, 'FaceColor', 'green');
end
hold on;
if exist('three', 'var')
    area(three.rainfall_excess.*1000, three.storage_runoff, 'FaceColor', 'blue');
end
if exist('thirty', 'var')
    area(thirty.rainfall_excess.*1000, thirty.storage_runoff,  'FaceColor', 'red');
end
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
if exist('ten', 'var')
    area(ten.rainfall_excess.*1000, 1000.*ten.storage_volume./(sum(~isnan(ten.dem(:)))*10*10), 'FaceColor', 'green');
end
hold on;
if exist('three', 'var')
    area(three.rainfall_excess.*1000, 1000.*three.storage_volume./(sum(~isnan(three.dem(:)))*3*3), 'FaceColor', 'blue');
end
if exist('thirty', 'var')
    area(thirty.rainfall_excess.*1000, 1000.*thirty.storage_volume./(sum(~isnan(thirty.dem(:)))*30*30),  'FaceColor', 'red');
end
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