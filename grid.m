close all
clear all


numberOfHearts   = 197;                        % number of hearts on a plate
volumePerHeart   = 5.4/numberOfHearts;        % volume per heart in mL

[Image_File, filepath] = uigetfile('*.jpg');           % import image
temperature      = 40;
volumetricFlowrate = 0.5 *2;
Grid_File        = imread(filepath);               % read file of coordinate grid
Coordinates      = uint32(csvread(Grid_File)); % predetermined coordinates of hearts
x                = Coordinates(:,1);
y                = Coordinates(:,2);
plate_size       = [1600 1060];                % height by width of plate

Image_Color      = imread(filepath);         
Plate_Color      = imcrop(Image_Color);
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
answer = questdlg('Are the circles where you would like to sample color?',...
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

    nondimensionalConcentration = zeros(1,numberOfHearts);
    subplot(3,3,3)
    for i = 1:numberOfHearts
        nondimensionalConcentration(1,i) = fsolve(@(x)0.085093.*x.^3 + ...
            -1.0526.*x.^2 + 1.7341 .*x+ 0.12132-heart_HSV(i,4),.5);
        clc
        plot(i,nondimensionalConcentration(1,i),'ko')
        hold on
    end
    xlabel('heart number')
    ylabel('C/C_{max}')
    
    count  = 0;
    for i = 1:length(nondimensionalConcentration)
        if nondimensionalConcentration(i) <= 0.05
            count = count+1;
        end
        if count >= 10
            heartOfCompletion =i;
            break
        end
    end
    
    volumeOfCompletion = volumePerHeart * heartOfCompletion;
    spaceTime = volumeOfCompletion / volumetricFlowrate;
    
    disp(['At a temperature of ',num2str(temperature),...
        ' C, flowrate of ', num2str(volumetricFlowrate),...
        ' mL/min the heart of 95% color reduction is ', num2str(heartOfCompletion), ...
        ' and the space time is ', num2str(spaceTime),' minutes.'])
    
elseif strcmp(answer,'Cancel')
    close all
    clear all
    clc
else
    grid
end