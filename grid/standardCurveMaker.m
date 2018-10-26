close all
clear all


%% Edit here

% Full filepath to save output to .txt file (CSV formatted)
dataSavingFilename = ['/Users/IsabelKaspriskie/Dropbox (MIT)/1026-S18'...
            '-T1-ContinuousFlow/Data/Experimental Data/Wittig Calibration2.txt'];


        

%% Do not change below this line.

numberOfHearts   = 197;                        % number of hearts on a plate
volumePerHeart   = 5.4/numberOfHearts;        % volume per heart in mL

% Have user browse for a file, from a specified "starting folder."
% For convenience in browsing, set a starting folder from which to browse.
startingFolder = ['/Users/IsabelKaspriskie/Dropbox (MIT)/'...
    '1026-S18-T1-ContinuousFlow/Data/Bales Camera/images-by'...
    '-category/Wittig calibration'];
if ~exist(startingFolder, 'dir')
	% If that folder doesn't exist, just start in the current folder.
	startingFolder = pwd;
end
% Get the name of the file that the user wants to use.
defaultFileName = fullfile(startingFolder, '*.*');
[baseFileName, folder] = uigetfile(defaultFileName, 'Select a file');
if baseFileName == 0
	% User clicked the Cancel button.
	return;
end
fullImageFilename = fullfile(folder, baseFileName);


dilution = str2double(extractBefore(baseFileName,'x'));
trial = extractBetween(baseFileName,'_','.JPG');

Grid_File        = 'coords.dat';               % read file of coordinate grid
Coordinates      = uint32(csvread(Grid_File)); % predetermined coordinates of hearts
x                = Coordinates(:,1);
y                = Coordinates(:,2);
plate_size       = [1600 1060];                % height by width of plate

Image_Color      = imread(fullImageFilename);         
Plate_Color      = imcrop(imrotate(Image_Color,90));
Plate_Sized      = imresize(Plate_Color, plate_size);
Plate_HSV        = rgb2hsv(Plate_Sized);

% Ask the user if the grid is mapped as they'd like it.
close all
checkGrid = figure('Position',[10,10,700,1000]);
imagesc(Plate_Sized);
hold on
for i = 1:numberOfHearts
    plot(x(i),y(i),'go')
    hold on
end
answer = questdlg('Continue?',...
    'Check grid');

if strcmp(answer, 'Yes')
    % Show the cropped image and where the grid thinks each heart center is
    close all
    
    figure('Position',[20,0,1400,600])
    subplot(3,3,[1,4,7])
    imagesc(Plate_Sized);
    hold on
    for i = 1:numberOfHearts
        plot(x(i),y(i),'go')
        hold on
    end

heart_HSV = zeros(197,5);
for i = 1:numberOfHearts
    
    heart_HSV(i,1) = x(i);
    heart_HSV(i,2) = y(i);
    
    % For these, the Plate_HSV takes the "height" pixel coordinate before
    % the "width," so must flip the y and x coordinates
    heart_HSV(i,3) = Plate_HSV(y(i),x(i),1); % Hue
    heart_HSV(i,4) = Plate_HSV(y(i),x(i),2); % Saturation
    heart_HSV(i,5) = Plate_HSV(y(i),x(i),3); % Value
end

subplot (3,3,2)
    for i = 1:numberOfHearts
        plot(i,heart_HSV(i,3),'ko') % Plot value as a function of heart number
        hold on
    end
    xlabel('heart number')
    ylabel('hue')
    subplot(3,3,5)
    for i = 1:numberOfHearts
        plot(i,heart_HSV(i,4),'ko') % Plot value as a function of heart number
        hold on
    end
    xlabel('heart number')
    ylabel('saturation')
    subplot(3,3,8)
    for i = 1:numberOfHearts
        plot(i,heart_HSV(i,5),'ko') % Plot value as a function of heart number
        hold on
    end
    xlabel('heart number')
    ylabel('value')
    
analysisOutput = [dilution, str2double(trial{1,1}), mean(heart_HSV(10:160,3)),...
    mean(heart_HSV(10:160,4))];


clc
disp(['dilution = ', num2str(analysisOutput(1)),'x'])
disp(['trial number = ', num2str(analysisOutput(2))])
disp(['average hue = ', num2str(analysisOutput(3))])
disp(['average saturation = ', num2str(analysisOutput(4))])

closingAnswer = questdlg('Continue?','Close plots');
    if strcmp(closingAnswer,'Yes')
        close all
            % Prompt user to save the data to a CSV file if they like the data.
        savingAnswer = questdlg(['Save to ' dataSavingFilename], ...
    'Save data');
        if strcmp(savingAnswer,'Yes')
            dlmwrite(dataSavingFilename,analysisOutput,'-append');
            disp('Saved.')
        else
            disp('Not saved.')
        end
    else

    end
    
elseif strcmp(answer,'Cancel')
    close all
    clear all
    clc
else
    grid
end