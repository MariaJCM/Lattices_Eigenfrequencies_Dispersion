%% Lattice Dispersion Visualization Toolkit
% Author: Maria Carrillo Munoz (Michigan Technological University)
% This script imports and visualizes eigenfrequency data from Abaqus .rpt files
% for a 2D rectangular lattice. It supports:
%   1. PATH analysis (boundary of IBZ, high-symmetry points)
%   2. GRID analysis (full IBZ interior, kx-ky grid)
% Both translational and rotational polarizations are visualized. Group velocity is plotted for specific types of propagation.
% 
% Input files: "example_data/abaqus_Frequencies_PATH.rpt" and "example_data/abaqus_Frequencies_GRID.rpt"
% (edit file paths below if needed)



%% Import the data
%% Initialize variables.
filename = '\abaqus_Frequencies_PATH.rpt';
delimiter = {'\t',' '};
startRow = 1;

%% Read columns of data as text:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
textscan(fileID, '%[^\n\r]', startRow-1, 'WhiteSpace', '', 'ReturnOnError', false, 'EndOfLine', '\r\n');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Convert the contents of columns containing numeric text to numbers.
% Replace non-numeric text with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
    % Converts text in the input cell array to numbers. Replaced non-numeric
    % text with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1)
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData(row), regexstr, 'names');
            numbers = result.numbers;

            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if numbers.contains(',')
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'))
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric text to numbers.
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


%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

%% Create output variable
FreqPATHRect_3_1 = cell2mat(raw);
ini=startRow+3;
%% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp R;

%% Data
[m1,m2] = size(FreqPATHRect_3_1);
NE=90; % number of eigenvalues (number of eigenfrequencies we calculated in abaqus).

i=ini:NE+3:m1;
j=2;      % frequency without any normalization
k=3;    % generalized mass
m_x=4; % mass Translation 
m_y=5; % mass Translation  
m_z=6; % mass Translation 
m_Ox=7; % mass rotation 
m_Oy=8; % mass rotation 
m_Oz=9; % mass rotation 
% 



n=30; % number of steps +1


for a=1:1:NE;

    F.(strcat('frequency',num2str(a)))=FreqPATHRect_3_1(i+a-1, j); % All the eigenfrequencies
    gmass.(strcat('frequency',num2str(a)))=FreqPATHRect_3_1(i+a-1, k);  % generalized mass
    em_x.(strcat('frequency',num2str(a)))=FreqPATHRect_3_1(i+a-1, m_x);  % effective mass x 
    em_y.(strcat('frequency',num2str(a)))=FreqPATHRect_3_1(i+a-1, m_y);  % effective mass y 
    em_z.(strcat('frequency',num2str(a)))=FreqPATHRect_3_1(i+a-1, m_z);  % effective mass z 
    em_Ox.(strcat('frequency',num2str(a)))=FreqPATHRect_3_1(i+a-1, m_Ox);  % effective mass Ox 
    em_Oy.(strcat('frequency',num2str(a)))=FreqPATHRect_3_1(i+a-1, m_Oy);  % effective mass Oy 
    em_Oz.(strcat('frequency',num2str(a)))=FreqPATHRect_3_1(i+a-1, m_Oz);  % effective mass Oz 



end

for a=1:2:(NE-1)
    etrans.(strcat('frequency',num2str(a)))=em_x.(strcat('frequency',num2str(a)))+em_x.(strcat('frequency',num2str(a+1)))+em_y.(strcat('frequency',num2str(a)))+em_y.(strcat('frequency',num2str(a+1)))+em_z.(strcat('frequency',num2str(a)))+em_z.(strcat('frequency',num2str(a+1)));
%     etransNORM.(strcat('frequency',num2str(a)))=etrans.(strcat('frequency',num2str(a)))/(3*1.6408609E-02);
    emass_x.(strcat('frequency',num2str(a)))=(em_x.(strcat('frequency',num2str(a)))+em_x.(strcat('frequency',num2str(a+1))))./etrans.(strcat('frequency',num2str(a)));
    emass_y.(strcat('frequency',num2str(a)))=(em_y.(strcat('frequency',num2str(a)))+em_y.(strcat('frequency',num2str(a+1))))./etrans.(strcat('frequency',num2str(a)));
    emass_z.(strcat('frequency',num2str(a)))=(em_z.(strcat('frequency',num2str(a)))+em_z.(strcat('frequency',num2str(a+1))))./etrans.(strcat('frequency',num2str(a)));

    erotac.(strcat('frequency',num2str(a)))=em_Ox.(strcat('frequency',num2str(a)))+em_Ox.(strcat('frequency',num2str(a+1)))+em_Oy.(strcat('frequency',num2str(a)))+em_Oy.(strcat('frequency',num2str(a+1)))+em_Oz.(strcat('frequency',num2str(a)))+em_Oz.(strcat('frequency',num2str(a+1)));
%     erotacNORM.(strcat('frequency',num2str(a)))=erotac.(strcat('frequency',num2str(a)))/(3*3.9225178E-07);
    emass_Ox.(strcat('frequency',num2str(a)))=(em_Ox.(strcat('frequency',num2str(a)))+em_Ox.(strcat('frequency',num2str(a+1))))./erotac.(strcat('frequency',num2str(a)));
    emass_Oy.(strcat('frequency',num2str(a)))=(em_Oy.(strcat('frequency',num2str(a)))+em_Oy.(strcat('frequency',num2str(a+1))))./erotac.(strcat('frequency',num2str(a)));
    emass_Oz.(strcat('frequency',num2str(a)))=(em_Oz.(strcat('frequency',num2str(a)))+em_Oz.(strcat('frequency',num2str(a+1))))./erotac.(strcat('frequency',num2str(a)));

end

% % % normale the frequencies if you want later to compare them with other lattices. 
% % FNor=F.frequency5(31); %norm by y motion
% % for a=1:1:NE;
% %     F.(strcat('frequency',num2str(a)))=F.(strcat('frequency',num2str(a)))/FNor;
% % end

total_ks=length(F.(strcat('frequency',num2str(a))));
k_vector=[0:1:total_ks-1]';

% % % % eff. mass accumulated if we want later to know how relevant our eigenmode is in the dynamic behavior of the lattice
% % % etransACU.(strcat('frequency',num2str(1)))=etrans.(strcat('frequency',num2str(1)));
% % % for a=3:2:NE; % to plot only the half of the eigenfrequencies.
% % % 
% % %     for b=1:1:length(emass_x.(strcat('frequency',num2str(a))));
% % %         etransACU.(strcat('frequency',num2str(a)))(b)=etransACU.(strcat('frequency',num2str(a-2)))(b)+etrans.(strcat('frequency',num2str(a)))(b);
% % %     end
% % % end

fig=0;




%% Frequency Dispersion Curves 

fig=fig+1;
figure(fig)
colormap (jet);


d=0.03;
for a=1:2:NE; % to plot only the half of the eigenfrequencies.
    plot(k_vector,F.(strcat('frequency',num2str(a))),'Color',[0.24705882370472 0.24705882370472 0.24705882370472]);
    hold on

%     for b=1+1:1:length(emass_x.(strcat('frequency',num2str(a)))); 
% 
% %        total mass in translation for each frequency , if we want to represent that in the dispersion curves
% 
%         scatter(k_vector(b),F.(strcat('frequency',num2str(a)))(b),30,etrans.(strcat('frequency',num2str(a)))(b)/(3*23.88057),'Marker','+','LineWidth',2.55); 
%         hold on

% %        total mass in translation ACUMULATED for each frequency, if we want to represent that information in the dispersion curves. 
%         
%         scatter(k_vector(b),F.(strcat('frequency',num2str(a)))(b),30,etransACU.(strcat('frequency',num2str(a)))(b)/(3*23.88057),'Marker','+','LineWidth',2.55); 
%         hold on

    % end     
end
% 
ylabel('Freq [Hz]');
% ylim([0.0 16]);
% yticks([0 2 4 6 8 10 12 14 16])
xlabel('k');
xlim([min(k_vector) max(k_vector)]);
xticks([min(k_vector)   total_ks/5   total_ks*2/5   total_ks*3/5   total_ks*4/5   max(k_vector)]);
xticklabels ({' O ',' X ',' M ',' O ', ' Y ', ' M '});
grid off
% grid minor
% colorbar
title('Frequency Dispersion Curves: Rectangular lattice 3/1, along the path OXMOYM');
ax = gca; % current axes
ax.FontSize = 22;
hold off


%% All the frequencies with mass_i factor 
fig=fig+1;
figure(fig)


for a=1:2:5%1:2:22%(NE-1); % to plot only some eigenfrequencies.
    % % plot(k_vector,F.(strcat('frequency',num2str(a))),'Color',[0.25 0.25 0.25]);
    % % hold on
    hd=0.5; % this is just to avoid overlapping of markers
    C=25; % this is a scale for all the markers; change according to the plot size to have a clear visualization

    for b=1:1:(length(emass_x.(strcat('frequency',num2str(a))))-1);  
        vd(b)=F.(strcat('frequency',num2str(a)))(b+1)-F.(strcat('frequency',num2str(a)))(b);
    end
    vd(length(emass_x.(strcat('frequency',num2str(a)))))=0;


    for b=1:1:length(emass_x.(strcat('frequency',num2str(a)))); 

        % s=C*(emass_x.(strcat('frequency',num2str(a)))(b)+0.0001);
        % scatter(k_vector(b),F.(strcat('frequency',num2str(a)))(b),s,emass_x.(strcat('frequency',num2str(a)))(b),'MarkerFaceAlpha',emass_x.(strcat('frequency',num2str(a)))(b),'MarkerFaceColor','blue','MarkerEdgeAlpha',emass_x.(strcat('frequency',num2str(a)))(b),'MarkerEdgeColor','blue','Marker','s','LineWidth',1.5) 
        % hold on
        % 
        % s=C*(emass_y.(strcat('frequency',num2str(a)))(b)+0.0001);
        % scatter(k_vector(b)+hd,F.(strcat('frequency',num2str(a)))(b)+0.5*vd(b), s ,emass_y.(strcat('frequency',num2str(a)))(b),'MarkerFaceAlpha',emass_y.(strcat('frequency',num2str(a)))(b),'MarkerFaceColor','green','MarkerEdgeAlpha',emass_y.(strcat('frequency',num2str(a)))(b),'MarkerEdgeColor','green','Marker','^','LineWidth',1.5) %
        % hold on
        % 
        % s=C*(emass_z.(strcat('frequency',num2str(a)))(b)+0.0001);
        % scatter(k_vector(b),F.(strcat('frequency',num2str(a)))(b), s ,emass_z.(strcat('frequency',num2str(a)))(b),'MarkerFaceAlpha',emass_z.(strcat('frequency',num2str(a)))(b),'MarkerFaceColor','red','MarkerEdgeAlpha',emass_z.(strcat('frequency',num2str(a)))(b),'MarkerEdgeColor','red','Marker','diamond','LineWidth',1.5) %
        % hold on
% 
        s=C*(emass_Ox.(strcat('frequency',num2str(a)))(b)+0.0001);
        scatter(k_vector(b),F.(strcat('frequency',num2str(a)))(b),s,emass_Ox.(strcat('frequency',num2str(a)))(b),'MarkerEdgeAlpha',emass_Ox.(strcat('frequency',num2str(a)))(b),'MarkerEdgeColor','blue','Marker','o','LineWidth',1.1) 
        hold on

        s=C*(emass_Oy.(strcat('frequency',num2str(a)))(b)+0.0001);
        scatter(k_vector(b)+hd,F.(strcat('frequency',num2str(a)))(b)+0.5*vd(b),s ,emass_Oy.(strcat('frequency',num2str(a)))(b),'MarkerEdgeAlpha',emass_Oy.(strcat('frequency',num2str(a)))(b),'MarkerEdgeColor','green','Marker','o','LineWidth',1.1) %
        hold on

        s=C*(emass_Oz.(strcat('frequency',num2str(a)))(b)+0.0001);
        scatter(k_vector(b),F.(strcat('frequency',num2str(a)))(b),s,emass_Oz.(strcat('frequency',num2str(a)))(b),'MarkerEdgeAlpha',emass_Oz.(strcat('frequency',num2str(a)))(b),'MarkerEdgeColor','red','Marker','o','LineWidth',1.1) %
        hold on
    end

end
% 
ylabel('\Omega');%('Freq [Hz]');  %Omega, or Freq in Hz, depending on having or not normalization.
% ylim([0.0 100]);
% yticks([0 2 4 6 8 10 12 14 16])
xlabel('k');
xlim([min(k_vector) max(k_vector)]);
xticks([min(k_vector)   total_ks/5   total_ks*2/5   total_ks*3/5   total_ks*4/5   max(k_vector)]);
xticklabels ({' O ',' X ',' M ',' O ', ' Y ', ' M '});
grid off%on
box on
% grid minor
% colorbar
title('Frequency Dispersion Curves- Effective Masses: Rectangular lattice 3/1, along the path OXMOYM 6dof');
ax = gca; % current axes
ax.FontSize = 22;
hold off





%% GRID

%% Import the data
%% Initialize variables.
filename = '...\abaqus_Frequencies_GRID.rpt';
delimiter = {'\t',' '};
startRow = 1;

%% Read columns of data as text:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
textscan(fileID, '%[^\n\r]', startRow-1, 'WhiteSpace', '', 'ReturnOnError', false, 'EndOfLine', '\r\n');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Convert the contents of columns containing numeric text to numbers.
% Replace non-numeric text with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
    % Converts text in the input cell array to numbers. Replaced non-numeric
    % text with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1)
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData(row), regexstr, 'names');
            numbers = result.numbers;

            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if numbers.contains(',')
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'))
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric text to numbers.
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


%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

%% Create output variable
FreqGRIDRect_3_1 = cell2mat(raw);
ini=startRow+3;
%% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp R;

%% Data
[m1,m2] = size(FreqGRIDRect_3_1);
NE=90; % number of eigenvalues

i=ini:NE+3:m1;
j=2;      % frequency without any normalization
k=3;    % generalized mass
m_x=4; % mass Translation 
m_y=5; % mass Translation  
m_z=6; % mass Translation 
m_Ox=7; % mass rotation 
m_Oy=8; % mass rotation 
m_Oz=9; % mass rotation 



n=31; % number of steps +1


FNor=1;
% % % normalizing the frequencies, if we want to compare then later with other lattices.
% FNor=F.frequency1(31);
% for a=1:1:NE;
%     F.(strcat('frequency',num2str(a)))=F.(strcat('frequency',num2str(a)))/FNor;
%  end


for a=1:1:NE;

    F.(strcat('frequency',num2str(a)))=FreqGRIDRect_3_1(i+a-1, j)/FNor; % All the eigenfrequencies
    gmass.(strcat('frequency',num2str(a)))=FreqGRIDRect_3_1(i+a-1, k);  % generalized mass
    em_x.(strcat('frequency',num2str(a)))=FreqGRIDRect_3_1(i+a-1, m_x);  % generalized mass x 
    em_y.(strcat('frequency',num2str(a)))=FreqGRIDRect_3_1(i+a-1, m_y);  % generalized mass y 
    em_z.(strcat('frequency',num2str(a)))=FreqGRIDRect_3_1(i+a-1, m_z);  % generalized mass z 
    em_Ox.(strcat('frequency',num2str(a)))=FreqGRIDRect_3_1(i+a-1, m_Ox);  % generalized mass Ox 
    em_Oy.(strcat('frequency',num2str(a)))=FreqGRIDRect_3_1(i+a-1, m_Oy);  % generalized mass Oy 
    em_Oz.(strcat('frequency',num2str(a)))=FreqGRIDRect_3_1(i+a-1, m_Oz);  % generalized mass Oz 

end

for a=1:2:(NE-1)
    etrans.(strcat('frequency',num2str(a)))=em_x.(strcat('frequency',num2str(a)))+em_x.(strcat('frequency',num2str(a+1)))+em_y.(strcat('frequency',num2str(a)))+em_y.(strcat('frequency',num2str(a+1)))+em_z.(strcat('frequency',num2str(a)))+em_z.(strcat('frequency',num2str(a+1)));

    emass_x.(strcat('frequency',num2str(a)))=(em_x.(strcat('frequency',num2str(a)))+em_x.(strcat('frequency',num2str(a+1))))./etrans.(strcat('frequency',num2str(a)));
    emass_y.(strcat('frequency',num2str(a)))=(em_y.(strcat('frequency',num2str(a)))+em_y.(strcat('frequency',num2str(a+1))))./etrans.(strcat('frequency',num2str(a)));
    emass_z.(strcat('frequency',num2str(a)))=(em_z.(strcat('frequency',num2str(a)))+em_z.(strcat('frequency',num2str(a+1))))./etrans.(strcat('frequency',num2str(a)));

    erotac.(strcat('frequency',num2str(a)))=em_Ox.(strcat('frequency',num2str(a)))+em_Ox.(strcat('frequency',num2str(a+1)))+em_Oy.(strcat('frequency',num2str(a)))+em_Oy.(strcat('frequency',num2str(a+1)))+em_Oz.(strcat('frequency',num2str(a)))+em_Oz.(strcat('frequency',num2str(a+1)));

    emass_Ox.(strcat('frequency',num2str(a)))=(em_Ox.(strcat('frequency',num2str(a)))+em_Ox.(strcat('frequency',num2str(a+1))))./erotac.(strcat('frequency',num2str(a)));
    emass_Oy.(strcat('frequency',num2str(a)))=(em_Oy.(strcat('frequency',num2str(a)))+em_Oy.(strcat('frequency',num2str(a+1))))./erotac.(strcat('frequency',num2str(a)));
    emass_Oz.(strcat('frequency',num2str(a)))=(em_Oz.(strcat('frequency',num2str(a)))+em_Oz.(strcat('frequency',num2str(a+1))))./erotac.(strcat('frequency',num2str(a)));

end





%%

Lx=0.5; % this is the dimension of this specific lattice
Ly=3.0*Lx; % This is the dimension of this lattice
Lz=0.0; % this is a 2D lattice
kx=0:pi/((n-1)*Lx):pi/Lx;
ky=0:pi/((n-1)*Ly):pi/Ly;
kz=0:pi/((n-1)*Lz):pi/Lz;
[kx_coef_MATRIX, ky_coef_MATRIX] = meshgrid( 0:pi/((n-1)*Lx):pi/Lx, 0:pi/((n-1)*Ly):pi/Ly );

for b=1:1:n;
    for c=1:1:n;

        for a=1:2:(NE-1);        
            F_MATRIX.(strcat('frequency',num2str(a)))(c,b)=F.(strcat('frequency',num2str(a)))(c+(b-1)*n);
            emass_x_MATRIX.(strcat('frequency',num2str(a)))(c,b)= emass_x.(strcat('frequency',num2str(a)))(c+(b-1)*n);
            emass_y_MATRIX.(strcat('frequency',num2str(a)))(c,b)= emass_y.(strcat('frequency',num2str(a)))(c+(b-1)*n);
            emass_z_MATRIX.(strcat('frequency',num2str(a)))(c,b)= emass_z.(strcat('frequency',num2str(a)))(c+(b-1)*n);
            emass_Ox_MATRIX.(strcat('frequency',num2str(a)))(c,b)= emass_Ox.(strcat('frequency',num2str(a)))(c+(b-1)*n);
            emass_Oy_MATRIX.(strcat('frequency',num2str(a)))(c,b)= emass_Oy.(strcat('frequency',num2str(a)))(c+(b-1)*n);
            emass_Oz_MATRIX.(strcat('frequency',num2str(a)))(c,b)= emass_Oz.(strcat('frequency',num2str(a)))(c+(b-1)*n);

        end

    end
end

fig=100;

%% All frequencies with mass_i factor  
%%% With the option alphaface I can make invisible the directions I do not
%%% want
for a=1:2:(NE-1);

    Ctr.(strcat('frequency',num2str(a)))(:,:,1) = emass_z_MATRIX.(strcat('frequency',num2str(a))); % red
    Ctr.(strcat('frequency',num2str(a)))(:,:,2) = emass_y_MATRIX.(strcat('frequency',num2str(a))); % green
    Ctr.(strcat('frequency',num2str(a)))(:,:,3) = emass_x_MATRIX.(strcat('frequency',num2str(a))); % blue

    Crot.(strcat('frequency',num2str(a)))(:,:,1) = emass_Oz_MATRIX.(strcat('frequency',num2str(a))); % red
    Crot.(strcat('frequency',num2str(a)))(:,:,2) = emass_Oy_MATRIX.(strcat('frequency',num2str(a))); % green
    Crot.(strcat('frequency',num2str(a)))(:,:,3) = emass_Ox_MATRIX.(strcat('frequency',num2str(a))); % blue
end

%% All directions together

fig=fig+1;
figure(fig)
for a=11:10:23%1:2:22%NE/2;
    transp_t=( (Ctr.(strcat('frequency',num2str(a)))(:,:,3)).^1+ (Ctr.(strcat('frequency',num2str(a)))(:,:,2).^1)+ (Ctr.(strcat('frequency',num2str(a)))(:,:,1)).^1  );
%     surf(kx_coef_MATRIX,ky_coef_MATRIX,F_MATRIX.(strcat('frequency',num2str(a))),Ctr.(strcat('frequency',num2str(a))),'FaceColor','flat','AlphaData', transp_t,'FaceAlpha','flat','EdgeColor',[0.5 0.5 0.5], 'EdgeAlpha','flat','LineStyle',':','LineWidth',0.1); %'EdgeAlpha',0.15,...
%     hold on 
    surf(kx_coef_MATRIX,ky_coef_MATRIX,F_MATRIX.(strcat('frequency',num2str(a))),Ctr.(strcat('frequency',num2str(a))),'FaceColor','flat','EdgeColor',[0.5 0.5 0.5], 'LineStyle',':','LineWidth',0.1); %'EdgeAlpha',0.15,...
    hold on 

%     transp_r=(  (Crot.(strcat('frequency',num2str(a)))(:,:,3)).^2 + (Crot.(strcat('frequency',num2str(a)))(:,:,2)).^2 + (Crot.(strcat('frequency',num2str(a)))(:,:,1)).^2  );
%     surf(kx_coef_MATRIX,ky_coef_MATRIX,F_MATRIX.(strcat('frequency',num2str(a))),Crot.(strcat('frequency',num2str(a))),'EdgeColor','flat','AlphaData',transp_r,'EdgeAlpha','flat','FaceColor','none'); %'EdgeAlpha',0.15,...
%     hold on 
%     mesh(kx_coef_MATRIX,ky_coef_MATRIX,F_MATRIX.(strcat('frequency',num2str(a))),Crot.(strcat('frequency',num2str(a))),'FaceLighting','gouraud','LineWidth',1.5,'FaceColor','none');%'FaceAlpha',0);
%     %     ,'MarkerEdgeColor',Crot,'MarkerSize',4,'Marker','o'
%     hold on
end
title('Rectangle 3/1 Lattice: effective mass_ i for any direction/rotation' )
% Label axes
xlabel kx
ylabel ky
zlabel Freq[Hz]
xlim([min(kx) max(kx)]); 
ylim([min(ky) max(ky)]);
% Create colorbar
% colorbar;
grid off
view( -30.7, 15.6 );
% view(-14.00, 0.30);
box on
ax = gca; % current axes
ax.FontSize = 22;
hold off

%% Specific directions
fig=fig+1;
figure(fig)
for a=1:2:36%NE;

%     surf(kx_coef_MATRIX,ky_coef_MATRIX,F_MATRIX.(strcat('frequency',num2str(a))),'AlphaData',Ctr.(strcat('frequency',num2str(a)))(:,:,3),'FaceAlpha','flat','FaceColor','blue', 'EdgeColor','none');
%     hold on
%     surf(kx_coef_MATRIX,ky_coef_MATRIX,F_MATRIX.(strcat('frequency',num2str(a))),'AlphaData',Ctr.(strcat('frequency',num2str(a)))(:,:,2),'FaceAlpha','flat','FaceColor','green', 'EdgeColor','none');
%     hold on
    surf(kx_coef_MATRIX,ky_coef_MATRIX,F_MATRIX.(strcat('frequency',num2str(a))),'AlphaData',Ctr.(strcat('frequency',num2str(a)))(:,:,1),'FaceAlpha','flat','FaceColor','red', 'EdgeColor','none');
    hold on


    surf(kx_coef_MATRIX,ky_coef_MATRIX,F_MATRIX.(strcat('frequency',num2str(a))),'FaceColor','none','LineWidth',0.15,'EdgeAlpha',0.55,'EdgeColor', [0.65 0.65 0.65]);
    hold on

%     surf(kx_coef_MATRIX,ky_coef_MATRIX,F_MATRIX.(strcat('frequency',num2str(a))),'AlphaData',Crot.(strcat('frequency',num2str(a)))(:,:,3),'EdgeAlpha','flat','EdgeColor','blue', 'FaceColor','none');
%     hold on
%     surf(kx_coef_MATRIX,ky_coef_MATRIX,F_MATRIX.(strcat('frequency',num2str(a))),'AlphaData',Crot.(strcat('frequency',num2str(a)))(:,:,2),'EdgeAlpha','flat','EdgeColor','green', 'FaceColor','none');
%     hold on
%     surf(kx_coef_MATRIX,ky_coef_MATRIX,F_MATRIX.(strcat('frequency',num2str(a))),'AlphaData',Crot.(strcat('frequency',num2str(a)))(:,:,1),'EdgeAlpha','flat','EdgeColor','red', 'FaceColor','none');
%     hold on
end
title('Rectangle 3 Lattice: effective mass_ i for any direction/rotation' )
% Label axes
xlabel kx
ylabel ky
zlabel Freq[Hz]
xlim([min(kx) max(kx)]);
ylim([min(ky) max(ky)]);
% Create colorbar
% colorbar;
grid off
view( -30.7, 15.6 );
box on
ax = gca; % current axes
ax.FontSize = 22;

hold off


%% CONSIDERING OX DIRECTION for the whole grid

%%
% 1. Gradient of each surface (independtly on the mode)

incry=ky_coef_MATRIX(2,1)-ky_coef_MATRIX(1,1);
incrx=kx_coef_MATRIX(1,2)-kx_coef_MATRIX(1,1);

for a=1:2:(NE-1);

    [Col,Row] = gradient(F_MATRIX.(strcat('frequency',num2str(a))));

    FX.(strcat('frequency',num2str(a)))=abs(Row/incrx);
    FY.(strcat('frequency',num2str(a)))=abs(Col/incry);
end


%% 2 Phi and Magnitude of each surface
for a=1:2:(NE-1);
    for b=1:1:length(kx);
        for c=1:1:length(ky);
            Phime.(strcat('frequency',num2str(a)))(b,c)=atan( FY.(strcat('frequency',num2str(a)))(b,c)/  FX.(strcat('frequency',num2str(a)))(b,c) );
            Magme.(strcat('frequency',num2str(a)))(b,c)= ( FY.(strcat('frequency',num2str(a)))(b,c)^2 + FX.(strcat('frequency',num2str(a)))(b,c)^2)^0.5;
        end
    end
end


%% 3 Phi, Magnitude and frequency of each surface for only longitudinal propagation in the X direction

for a=1:2:(NE-1);
    Phime_LONGX.(strcat('frequency',num2str(a)))=zeros(n,n);
    F_MATRIX_LONGX.(strcat('frequency',num2str(a)))=zeros(n,n);
    Magme_LONGX.(strcat('frequency',num2str(a)))=zeros(n,n);
    for b=1:1:length(kx);
        for c=1:1:length(ky);

            if emass_x_MATRIX.(strcat('frequency',num2str(a)))(b,c)>=0.9 && emass_y_MATRIX.(strcat('frequency',num2str(a)))(b,c)<=0.15 && emass_z_MATRIX.(strcat('frequency',num2str(a)))(b,c)<=0.1 ;
                if F_MATRIX.(strcat('frequency',num2str(a)))(b,c)<=16; %this is to represent only values under this one
                    Phime_LONGX.(strcat('frequency',num2str(a)))(b,c)=Phime.(strcat('frequency',num2str(a)))(b,c);
                    F_MATRIX_LONGX.(strcat('frequency',num2str(a)))(b,c)=F_MATRIX.(strcat('frequency',num2str(a)))(b,c);
                    Magme_LONGX.(strcat('frequency',num2str(a)))(b,c)=Magme.(strcat('frequency',num2str(a)))(b,c);
                end
            end  

        end

    end

end

%%
%%%% Plot of the Group Velocity for longitudinal propagation in the X  
%%%% direction

fig=fig+1;
figure(fig);
colormap (jet);

for a=1:2:(NE-1);

    for b=1:1:length(kx);
        for c=1:1:length(ky);

                polarscatter( Phime_LONGX.(strcat('frequency',num2str(a)))(b,c), F_MATRIX_LONGX.(strcat('frequency',num2str(a)))(b,c),500, Magme_LONGX.(strcat('frequency',num2str(a)))(b,c), '.'   );
        end  
        hold on           
    end
    hold on
end


hold off
thetalim([0 90])
ax = gca; % current axes

ax.FontSize = 18;%22
ax.GridAlpha = 0.5;
ax.MinorGridAlpha = 0.5;
title('\fontsize{14}RECTANGULAR Lattice : Group Velocity of longitudinal propagation in the X direction' )
% Create colorbar
colorbar;



%% 4 Phi, Magnitude and frequency of each surface for only shear in plane respect to the X direction (so it is y displacement)

for a=1:2:(NE-1);
    Phime_SHEARX.(strcat('frequency',num2str(a)))=zeros(n,n);
    F_MATRIX_SHEARX.(strcat('frequency',num2str(a)))=zeros(n,n);
    Magme_SHEARX.(strcat('frequency',num2str(a)))=zeros(n,n);
    for b=1:1:length(kx);
        for c=1:1:length(ky);

            if emass_y_MATRIX.(strcat('frequency',num2str(a)))(b,c)>=0.9 && emass_x_MATRIX.(strcat('frequency',num2str(a)))(b,c)<=0.15 && emass_z_MATRIX.(strcat('frequency',num2str(a)))(b,c)<=0.1 ;
                if F_MATRIX.(strcat('frequency',num2str(a)))(b,c)<=16; %this is to represent only values under this one
                    Phime_SHEARX.(strcat('frequency',num2str(a)))(b,c)=Phime.(strcat('frequency',num2str(a)))(b,c);
                    F_MATRIX_SHEARX.(strcat('frequency',num2str(a)))(b,c)=F_MATRIX.(strcat('frequency',num2str(a)))(b,c);
                    Magme_SHEARX.(strcat('frequency',num2str(a)))(b,c)=Magme.(strcat('frequency',num2str(a)))(b,c);
                end
            end  

        end

    end

end

%%
%%%% Plot of the Group Velocity for shear in plane respect to the X  
%%%% direction

fig=fig+1;
figure(fig);
colormap (jet);

for a=1:2:(NE-1);

    for b=1:1:length(kx);
        for c=1:1:length(ky);

                polarscatter( Phime_SHEARX.(strcat('frequency',num2str(a)))(b,c), F_MATRIX_SHEARX.(strcat('frequency',num2str(a)))(b,c),500, Magme_SHEARX.(strcat('frequency',num2str(a)))(b,c), '.'   );
        end  
        hold on           
    end
    hold on
end


hold off
thetalim([0 90])
ax = gca; % current axes

ax.FontSize = 18;%22
ax.GridAlpha = 0.5;
ax.MinorGridAlpha = 0.5;
title('\fontsize{14}SQUARE Lattice : Group Velocity of shear in plane propagation respect to the X direction' )
% Create colorbar
colorbar;





%% 5 Phi, Magnitude and frequency of each surface for only out of plane 

for a=1:2:(NE-1);
    Phime_OUT.(strcat('frequency',num2str(a)))=zeros(n,n);
    F_MATRIX_OUT.(strcat('frequency',num2str(a)))=zeros(n,n);
    Magme_OUT.(strcat('frequency',num2str(a)))=zeros(n,n);
    for b=1:1:length(kx);
        for c=1:1:length(ky);

            if emass_z_MATRIX.(strcat('frequency',num2str(a)))(b,c)>=0.35 && emass_x_MATRIX.(strcat('frequency',num2str(a)))(b,c)<=0.15 && emass_y_MATRIX.(strcat('frequency',num2str(a)))(b,c)<=0.15 ;
                if F_MATRIX.(strcat('frequency',num2str(a)))(b,c)<=16; %this is to represent only values under this one
                    Phime_OUT.(strcat('frequency',num2str(a)))(b,c)=Phime.(strcat('frequency',num2str(a)))(b,c);
                    F_MATRIX_OUT.(strcat('frequency',num2str(a)))(b,c)=F_MATRIX.(strcat('frequency',num2str(a)))(b,c);
                    Magme_OUT.(strcat('frequency',num2str(a)))(b,c)=Magme.(strcat('frequency',num2str(a)))(b,c);
                end
            end  

        end

    end

end

%%
%%%% Plot of the Group Velocity for shear out of plane  
%%%% direction

fig=fig+1;
figure(fig);
colormap (jet);

for a=1:2:(NE-1);

    for b=1:1:length(kx);
        for c=1:1:length(ky);

                polarscatter( Phime_OUT.(strcat('frequency',num2str(a)))(b,c), F_MATRIX_OUT.(strcat('frequency',num2str(a)))(b,c),500, Magme_OUT.(strcat('frequency',num2str(a)))(b,c), '.'   );
        end  
        hold on           
    end
    hold on
end


hold off
thetalim([0 90])
ax = gca; % current axes

ax.FontSize = 18;%22
ax.GridAlpha = 0.5;
ax.MinorGridAlpha = 0.5;
title('\fontsize{14}SQUARE Lattice : Group Velocity of out of plane propagation' )
% Create colorbar
colorbar;
