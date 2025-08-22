%% Rectangular Lattice Dispersion Visualization Toolkit
% Author: Maria Carrillo Munoz (Michigan Technological University)
% This script imports and visualizes eigenfrequency data from Abaqus .rpt files
% for a 2D rectangular lattice. It supports:
%   1. PATH analysis (boundary of IBZ, high-symmetry points)
%   2. GRID analysis (full IBZ interior, kx-ky grid)
% Both translational and rotational polarizations are visualized.
% 
% Input files: "example_data/abaqus_Frequencies_PATH.rpt" and "example_data/abaqus_Frequencies_GRID.rpt"
% (edit file paths below if needed)

%% ===========================================================
%  PART 1: PATH ANALYSIS (BOUNDARY OF IBZ)
% ===========================================================

% ---- 1.1 Import PATH Data ----
filename_path = 'example_data/abaqus_Frequencies_PATH.rpt';
delimiter = {'\t',' '};
startRow = 1;
formatSpec = repmat('%s',1,16) + "%[^\n\r]"; % 16 columns

fileID = fopen(filename_path,'r');
textscan(fileID, '%[^\n\r]', startRow-1, 'WhiteSpace', '', 'ReturnOnError', false, 'EndOfLine', '\r\n');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'ReturnOnError', false);
fclose(fileID);

% ---- 1.2 Clean and Convert Data ----
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));
for col=1:16
    rawData = dataArray{col};
    for row=1:size(rawData, 1)
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData(row), regexstr, 'names');
            numbers = result.numbers;
            invalidThousandsSeparator = false;
            if contains(numbers, ',')
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'))
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            if ~invalidThousandsSeparator
                numbers = textscan(char(strrep(numbers, ',', '')), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch
            raw{row, col} = rawData{row};
        end
    end
end
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells
FreqPATH = cell2mat(raw); % Final numeric matrix
clearvars filename_path delimiter startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp R;

% ---- 1.3 Structure Data for Analysis ----
[m1,~] = size(FreqPATH);
NE = 90; % Number of eigenvalues (modes) per k-point
ini = startRow + 3;
i = ini:NE+3:m1; % Start of each block

% Indices for physical quantities (edit if column order changes!)
j_freq = 2; m_x=4; m_y=5; m_z=6; m_Ox=7; m_Oy=8; m_Oz=9;
n_path = 30; % Number of k-points along path

for a = 1:NE
    F_path(:,a) = FreqPATH(i+a-1, j_freq);        % Frequency (Hz)
    em_x_path(:,a) = FreqPATH(i+a-1, m_x);        % Effective mass X
    em_y_path(:,a) = FreqPATH(i+a-1, m_y);        % Effective mass Y
    em_z_path(:,a) = FreqPATH(i+a-1, m_z);        % Effective mass Z
    em_Ox_path(:,a) = FreqPATH(i+a-1, m_Ox);      % Rotational mass Ox
    em_Oy_path(:,a) = FreqPATH(i+a-1, m_Oy);      % Rotational mass Oy
    em_Oz_path(:,a) = FreqPATH(i+a-1, m_Oz);      % Rotational mass Oz
end

% ---- 1.4 Normalize Effective Masses ----
% Normalization reveals the dominant vibration direction for each mode
for a=1:2:(NE-1)
    etrans = em_x_path(:,a) + em_x_path(:,a+1) + em_y_path(:,a) + em_y_path(:,a+1) + em_z_path(:,a) + em_z_path(:,a+1);
    emass_x_path(:,a) = (em_x_path(:,a) + em_x_path(:,a+1)) ./ etrans;
    emass_y_path(:,a) = (em_y_path(:,a) + em_y_path(:,a+1)) ./ etrans;
    emass_z_path(:,a) = (em_z_path(:,a) + em_z_path(:,a+1)) ./ etrans;
    erotac = em_Ox_path(:,a) + em_Ox_path(:,a+1) + em_Oy_path(:,a) + em_Oy_path(:,a+1) + em_Oz_path(:,a) + em_Oz_path(:,a+1);
    emass_Ox_path(:,a) = (em_Ox_path(:,a) + em_Ox_path(:,a+1)) ./ erotac;
    emass_Oy_path(:,a) = (em_Oy_path(:,a) + em_Oy_path(:,a+1)) ./ erotac;
    emass_Oz_path(:,a) = (em_Oz_path(:,a) + em_Oz_path(:,a+1)) ./ erotac;
end

k_vector = 0:(n_path-1);

% ---- 1.5 Plot: Dispersion Curves (Translations) ----
figure; hold on;
for a=1:2:NE
    plot(k_vector, F_path(:,a), 'Color', [0.2 0.2 0.2]);
end
xlabel('Wave vector (path index)');
ylabel('Frequency [Hz]');
title('Dispersion Curves - Translations (Path)');
xticks([0 n_path/5 n_path*2/5 n_path*3/5 n_path*4/5 n_path-1]);
xticklabels({'O','X','M','O','Y','M'});
set(gca,'FontSize',14); box on; hold off;

% ---- 1.6 Plot: Dispersion Curves (Rotations) ----
figure; hold on;
for a=1:2:NE
    plot(k_vector, F_path(:,a), 'Color', [0.6 0.1 0.8]);
end
xlabel('Wave vector (path index)');
ylabel('Frequency [Hz]');
title('Dispersion Curves - Rotations (Path)');
xticks([0 n_path/5 n_path*2/5 n_path*3/5 n_path*4/5 n_path-1]);
xticklabels({'O','X','M','O','Y','M'});
set(gca,'FontSize',14); box on; hold off;

% ---- 1.7 Plot: Polarization Markers (Translations & Rotations) ----
figure; hold on;
colors = {'b','g','r'}; % X, Y, Z
rot_colors = {'c','m','y'}; % Ox, Oy, Oz
for a=1:2:5
    for b=1:n_path
        % Translation markers
        scatter(k_vector(b), F_path(b,a), 40*emass_x_path(b,a), colors{1}, 'filled', 'Marker', 's');
        scatter(k_vector(b)+0.5, F_path(b,a), 40*emass_y_path(b,a), colors{2}, 'filled', 'Marker', '^');
        scatter(k_vector(b), F_path(b,a), 40*emass_z_path(b,a), colors{3}, 'filled', 'Marker', 'd');
        % Rotation markers
        scatter(k_vector(b), F_path(b,a), 40*emass_Ox_path(b,a), rot_colors{1}, 'filled', 'Marker', 'o');
        scatter(k_vector(b)+0.5, F_path(b,a), 40*emass_Oy_path(b,a), rot_colors{2}, 'filled', 'Marker', 'o');
        scatter(k_vector(b), F_path(b,a), 40*emass_Oz_path(b,a), rot_colors{3}, 'filled', 'Marker', 'o');
    end
end
xlabel('Wave vector (path index)');
ylabel('Frequency [Hz]');
title('Polarization Markers (Translations & Rotations)');
xticks([0 n_path/5 n_path*2/5 n_path*3/5 n_path*4/5 n_path-1]);
xticklabels({'O','X','M','O','Y','M'});
set(gca,'FontSize',14); box on; hold off;

%% ===========================================================
%  PART 2: GRID ANALYSIS (FULL IBZ INTERIOR)
% ===========================================================

% ---- 2.1 Import GRID Data ----
filename_grid = 'example_data/abaqus_Frequencies_GRID.rpt';
fileID = fopen(filename_grid,'r');
textscan(fileID, '%[^\n\r]', startRow-1, 'WhiteSpace', '', 'ReturnOnError', false, 'EndOfLine', '\r\n');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'ReturnOnError', false);
fclose(fileID);

raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));
for col=1:16
    rawData = dataArray{col};
    for row=1:size(rawData, 1)
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData(row), regexstr, 'names');
            numbers = result.numbers;
            invalidThousandsSeparator = false;
            if contains(numbers, ',')
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'))
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            if ~invalidThousandsSeparator
                numbers = textscan(char(strrep(numbers, ',', '')), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch
            raw{row, col} = rawData{row};
        end
    end
end
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw);
raw(R) = {NaN};
FreqGRID = cell2mat(raw);
clearvars filename_grid fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp R;

% ---- 2.2 Structure Data for Analysis ----
[m1,~] = size(FreqGRID);
NE = 90; % Number of eigenvalues
ini = startRow + 3;
i = ini:NE+3:m1;
n_grid = 31; % kx, ky points per axis

for a = 1:NE
    F_grid(:,a) = FreqGRID(i+a-1, j_freq);
    em_x_grid(:,a) = FreqGRID(i+a-1, m_x);
    em_y_grid(:,a) = FreqGRID(i+a-1, m_y);
    em_z_grid(:,a) = FreqGRID(i+a-1, m_z);
    em_Ox_grid(:,a) = FreqGRID(i+a-1, m_Ox);
    em_Oy_grid(:,a) = FreqGRID(i+a-1, m_Oy);
    em_Oz_grid(:,a) = FreqGRID(i+a-1, m_Oz);
end

% ---- 2.3 Normalize Effective Masses (GRID) ----
for a=1:2:(NE-1)
    etrans = em_x_grid(:,a) + em_x_grid(:,a+1) + em_y_grid(:,a) + em_y_grid(:,a+1) + em_z_grid(:,a) + em_z_grid(:,a+1);
    emass_x_grid(:,a) = (em_x_grid(:,a) + em_x_grid(:,a+1)) ./ etrans;
    emass_y_grid(:,a) = (em_y_grid(:,a) + em_y_grid(:,a+1)) ./ etrans;
    emass_z_grid(:,a) = (em_z_grid(:,a) + em_z_grid(:,a+1)) ./ etrans;
    erotac = em_Ox_grid(:,a) + em_Ox_grid(:,a+1) + em_Oy_grid(:,a) + em_Oy_grid(:,a+1) + em_Oz_grid(:,a) + em_Oz_grid(:,a+1);
    emass_Ox_grid(:,a) = (em_Ox_grid(:,a) + em_Ox_grid(:,a+1)) ./ erotac;
    emass_Oy_grid(:,a) = (em_Oy_grid(:,a) + em_Oy_grid(:,a+1)) ./ erotac;
    emass_Oz_grid(:,a) = (em_Oz_grid(:,a) + em_Oz_grid(:,a+1)) ./ erotac;
end

% ---- 2.4 Reshape Data to [n_grid x n_grid x modes] ----
F_grid_reshaped = zeros(n_grid,n_grid,NE);
emass_x_grid_reshaped = zeros(n_grid,n_grid,NE);
emass_y_grid_reshaped = zeros(n_grid,n_grid,NE);
emass_z_grid_reshaped = zeros(n_grid,n_grid,NE);
emass_Ox_grid_reshaped = zeros(n_grid,n_grid,NE);
emass_Oy_grid_reshaped = zeros(n_grid,n_grid,NE);
emass_Oz_grid_reshaped = zeros(n_grid,n_grid,NE);

for a=1:2:(NE-1)
    F_grid_reshaped(:,:,a) = reshape(F_grid(:,a), n_grid, n_grid);
    emass_x_grid_reshaped(:,:,a) = reshape(emass_x_grid(:,a), n_grid, n_grid);
    emass_y_grid_reshaped(:,:,a) = reshape(emass_y_grid(:,a), n_grid, n_grid);
    emass_z_grid_reshaped(:,:,a) = reshape(emass_z_grid(:,a), n_grid, n_grid);
    emass_Ox_grid_reshaped(:,:,a) = reshape(emass_Ox_grid(:,a), n_grid, n_grid);
    emass_Oy_grid_reshaped(:,:,a) = reshape(emass_Oy_grid(:,a), n_grid, n_grid);
    emass_Oz_grid_reshaped(:,:,a) = reshape(emass_Oz_grid(:,a), n_grid, n_grid);
end

% ---- 2.5 Set up kx, ky grid ----
Lx = 0.5; Ly = 3.0*Lx; % Unit cell sizes
kx = linspace(0, pi/Lx, n_grid);
ky = linspace(0, pi/Ly, n_grid);
[kx_grid_mesh, ky_grid_mesh] = meshgrid(kx, ky);

% ---- 2.6 Plot: Dispersion Surfaces (Translations) ----
figure; hold on;
for a=11:10:23
    surf(kx_grid_mesh, ky_grid_mesh, F_grid_reshaped(:,:,a), cat(3, emass_z_grid_reshaped(:,:,a), emass_y_grid_reshaped(:,:,a), emass_x_grid_reshaped(:,:,a)), ...
        'FaceColor', 'interp', 'EdgeColor', 'none');
end
xlabel('k_x'); ylabel('k_y'); zlabel('Frequency [Hz]');
title('Dispersion Surfaces - Translations (Grid)');
set(gca,'FontSize',14); view(-30,18); colorbar; box on; hold off;

% ---- 2.7 Plot: Dispersion Surfaces (Rotations) ----
figure; hold on;
for a=11:10:23
    surf(kx_grid_mesh, ky_grid_mesh, F_grid_reshaped(:,:,a), cat(3, emass_Oz_grid_reshaped(:,:,a), emass_Oy_grid_reshaped(:,:,a), emass_Ox_grid_reshaped(:,:,a)), ...
        'FaceColor', 'interp', 'EdgeColor', 'none');
end
xlabel('k_x'); ylabel('k_y'); zlabel('Frequency [Hz]');
title('Dispersion Surfaces - Rotations (Grid)');
set(gca,'FontSize',14); view(-30,18); colorbar; box on; hold off;

disp('Analysis complete. See figures for results.');