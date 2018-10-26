% Isabel Kaspriskie (ikaspriskie@gmail.com)
% Team oneFlo: 10.26 Chemical Engineering Projects Lab, Spring 2018
% Image analysis code for Corning Advanced-Flow Reactor (AFR)


close all
clear all
%% Settings

reaction = "Wittig"; % Wittig or bromination
dataSavingFilename = ['/Users/IsabelKaspriskie/Dropbox (MIT)/1026-S18'...
            '-T1-ContinuousFlow/Data/Experimental Data/Wittig Kinetics.txt'];
                          % where to output data for further processing

                          
                          
                          
                          
%% Do not change below this line

numberOfHearts = 197;     % number of hearts on a plate

% Have user browse for a file, from a specified "starting folder."
% For convenience in browsing, set a starting folder from which to browse.
startingFolder = ['/Users/IsabelKaspriskie/Dropbox (MIT)/'...
    '1026-S18-T1-ContinuousFlow/Data/Bales Camera/images-by'...
    '-category'];
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

volumePerHeart          = 5.4/numberOfHearts;     % volume per heart in mL

tempGuess               = extractBefore(baseFileName,'_');
flowrateGuess           = extractBetween(baseFileName,'_','_');
numberGuess             = extractBefore(extractAfter(...
    extractAfter(baseFileName, '_'),'_'),'.J|.j');

% Check that the input temperature, flowrate, and picture number are
% correct.
prompt                  = {'Temperature [Celcius]:',...
                            'Single pump flowrate [mL/min]:',...
                            'Picture Number: '};
title                   = baseFileName;
dims                    = [1 35];
definput                = {tempGuess,flowrateGuess{1,1},numberGuess};
answer                  = inputdlg(prompt,title,dims,definput);

temperature             = str2double(answer{1,1});
totalVolumetricFlowrate = str2double(answer{2,1}) *2;
trialNumber             = str2double(answer{3,1});
Grid_File               = 'coords.dat';               % read file of coordinate grid
Coordinates             = uint32(csvread(Grid_File)); % predetermined coordinates of hearts
x                       = Coordinates(:,1);
y                       = Coordinates(:,2);
plate_size              = [1600 1060];                % height by width of plate

Image_Color             = imread(fullImageFilename);
Plate_Color             = imcrop(imrotate(Image_Color,90));
Plate_Sized             = imresize(Plate_Color, plate_size);
Plate_HSV               = rgb2hsv(Plate_Sized);
Plate_Gray              = rgb2gray(Plate_Sized);


% Ask the user if the grid is mapped as they'd like it.
close all
checkGrid = figure('Position',[10,10,650,1000]);
imshow(Plate_Sized);
hold on
for i = 1:numberOfHearts
    plot(x(i),y(i),'go')
    hold on
end
checkGridAnswer = questdlg('Are the circles where you would like to sample color?',...
    'Check grid');

if strcmp(checkGridAnswer, 'Yes')
    % Show the cropped image and where the grid thinks each heart center is
    close all
    
    figure('Position',[20,0,1400,600])
    subplot(3,3,[1,4,7])
    imshow(Plate_Sized);
    hold on
    for i = 1:numberOfHearts
        plot(x(i),y(i),'go')
        hold on
    end
    xlabel([num2str(temperature) 'ºC, ' num2str(totalVolumetricFlowrate) 'mL/min' ])

    heart_values = zeros(197,6);
    for i = 1:numberOfHearts
        heart_values(i,1) = x(i);
        heart_values(i,2) = y(i);

        % For these, the Plate_HSV takes the "height" pixel coordinate before
        % the "width," so must flip the y and x coordinates
        heart_values(i,3) = Plate_HSV(y(i),x(i),1); % Hue
        heart_values(i,4) = Plate_HSV(y(i),x(i),2); % Saturation
        heart_values(i,5) = Plate_HSV(y(i),x(i),3); % Value
        heart_values(i,6) = Plate_Gray(y(i),x(i)); % Gray
    end

    subplot (3,3,2)
    for i = 1:numberOfHearts
        plot(i,heart_values(i,3),'ko') % Plot hue as a function of heart number
        hold on
    end
    xlabel('heart number')
    ylabel('hue')
    subplot(3,3,5)
    for i = 1:numberOfHearts
        plot(i,heart_values(i,4),'ko') % Plot saturation as a function of heart number
        hold on
    end
    xlabel('heart number')
    ylabel('saturation')
    subplot(3,3,8)
    for i = 1:numberOfHearts
        plot(i,heart_values(i,5),'ko') % Plot value as a function of heart number
        hold on
    end
    xlabel('heart number')
    ylabel('value')

    nondimensionalConcentration = zeros(1,numberOfHearts);
    time = zeros(numberOfHearts,1);
    subplot(3,3,3)
    for i = 1:numberOfHearts
        if strcmp(reaction, "Wittig")
            heart = i;
            nondimensionalConcentration(1,i) = fsolve(@(x)0.085093.*x.^3 + ...
                -1.0526.*x.^2 + 1.7341 .*x+ 0.12132-heart_values(i,4),.5);
            clc
            time(i) = volumePerHeart * heart / totalVolumetricFlowrate*60; 
            plot(time(i),nondimensionalConcentration(1,i),'ko');
            hold on
        elseif strcmp(reaction,"bromination")
            heart = i
            nondimensionalConcentration(1,i) = fsolve(@(x) -0.2189.*x.^3 + ...
                0.5056.*x.^2-0.4073.*x+0.1505-heart_values(i,3),.5);
            clc
            time(i) = volumePerHeart * heart / totalVolumetricFlowrate*60; 
            plot(time(i),nondimensionalConcentration(1,i),'ko');
            hold on
        end
    end
    xlabel('time [s]')
    ylabel('C/C_o')
    
f = fit(time,nondimensionalConcentration','exp1');
coeffvals = coeffvalues(f);
rateConstant = -coeffvals(2);
subplot(3,3,6)
fplot(@(x)exp(-rateConstant.*x),[0 350],'linewidth',2,'color','k');
xlabel('time [s]')
ylabel('C/C_o')
legend(['k = ' num2str(round(rateConstant,3,'significant')) ' s^{-1}'])

subplot(3,3,9)
plot(f,time,nondimensionalConcentration')
xlabel('time [s]')
ylabel('C/C_o')

    
%     subplot(3,3,6)
%     for i = 1:numberOfHearts
%         time(i) = volumePerHeart * i / totalVolumetricFlowrate;
%         plot(time(i),heart_values(i,6),'ko') % Plot value as a function of heart number
%         hold on
%     end
%     axis([0 1.5 0 200])
%     xlabel('time [min]')
%     ylabel('grayscale value')
%     
    count  = 0;
    for i = 1:length(nondimensionalConcentration)
        if nondimensionalConcentration(i) <= 0.05
            count = count+1;
        end
        if count >= 10
            heartOfCompletion =i;
            break
        else
            heartOfCompletion = nan;
        end
    end
    
    volumeOfCompletion = heartOfCompletion * volumePerHeart;
    timeOfCompletion = volumeOfCompletion/totalVolumetricFlowrate; % mL / (mL/min)
   
    
%     Matrix to output to CSV. 
%     Trial number, temperature [C], single pump
%     flowrate [mL/min], heart of completion, volume of completion [mL],
%     coefficient of exponential fit, rate constant [s^-1].

    analysisOutput = [trialNumber, temperature, totalVolumetricFlowrate,...
        volumeOfCompletion,heartOfCompletion,rateConstant];
    
    disp(['trial: ',num2str(analysisOutput(1))])
    disp(['temperature [K]: ',num2str(analysisOutput(2))])
    disp(['volumetric flowrate [mL/min]: ',num2str(analysisOutput(3))])
    disp(['rate constant [1/s]: ',num2str(analysisOutput(6))])

    closingAnswer = questdlg('Continue?','Close plots');
    if strcmp(closingAnswer,'Yes')
        close all
            % Prompt user to save the data to a CSV file if they like the data.
        savingAnswer = questdlg('Would you like to save this information for post-processing?', ...
    'Save data');
        if strcmp(savingAnswer,'Yes')
            dlmwrite(dataSavingFilename,analysisOutput,'-append');
            disp('Saved.')
        else
            disp('Cancelled.')
        end
    else

    end

       
elseif strcmp(checkGridAnswer,'Cancel')
    close all
    clear all
    clc
else
    grid
end