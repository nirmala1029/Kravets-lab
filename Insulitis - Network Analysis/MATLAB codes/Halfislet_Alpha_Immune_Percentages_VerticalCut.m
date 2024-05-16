%% This code slices an islet into two halves using a vertical line, and computes the percentage of islet cells and immune cells in each half
%% Credit - Nirmala Balasenthilkumaran, Vira Kravets
%% Last modified - April 2024

%% Clear command window, and all saved variables

clc
clear all 
close all

%% Initialization  

filepath = 'C:\Nirmala\UCSD\Kravets Lab\T1D Project\All islets\'; %Initialize filepath to load files containing islet masks and the positions of various cell types in the islet
sno = 0; %Initialize looping variable
scaling = 0.4964671; %Relationship between 1 pixel and 1 um

%% Percentage analysis of user-selected islets

for isletno = 1:134 %Update this line with a list of islets to be analyzed

    sno = sno + 1; %Update looping variable

    mask = logical(imread([filepath,['islet ', num2str(isletno), ' mask 20.tif']])); %Load the mask encompassing islets cells and immune cells, 20 um from the islet rim

    %% Create a bounding box for the mask
    
    % Initialize variables used to store the row and column numbers of the edges of the mask

    uppervert = 0;
    lowervert = 0;
    leftvert = 0;
    rightvert = 0;

    % Determine the uppermost point of the mask

    for i = 1:size(mask,1) %Iterate through the different rows of the mask, from the first row to the last row
        if checkforones(mask(i,:)) %The index of the first non-zero pixel containing row is stored
            uppervert = i;
            break;
        end
    end

    % Determine the lowermost point of the mask

    for i = size(mask,1):-1:1 %Iterate through the different rows of the mask, from the last row to the first row
        if checkforones(mask(i,:)) %The index of the first non-zero pixel containing row is stored
            lowervert = i;
            break;
        end
    end

    % Determine the leftmsot point of the mask

    for i = 1:size(mask,2) %Iterate through the different columns of the mask, from the first column to the last column
        if checkforones(mask(:,i)) %The index of the first non-zero pixel containing column is stored
            leftvert = i;
            break;
        end
    end

    % Determine the right most point of the mask

    for i = size(mask,2):-1:1 %Iterate through the different columns of the mask, from the last column to the first column
        if checkforones(mask(:,i)) %The index of the first non-zero pixel containing column is stored
            rightvert = i;
            break;
        end
    end

    % Generate a mask using the edges of the bounding box.

    bbox = zeros(size(mask));
    bbox(uppervert:lowervert,leftvert:rightvert) = 1;

    %% Split the mask into two halves

    % Determine the dimensions of the mask

    width = rightvert - leftvert;
    height = lowervert - uppervert; 

    % Determine the half length and half width of the mask

    partwidth = floor(width/2);
    partheight = floor(height/2);

    % Initialize cells two store the two halves of the mask

    bboxmask = cell(1,2);
    halfmask = cell(1,2);

    for i=1:2
        halfmask{i} = zeros(size(mask));
        bboxmask{i} = zeros(size(mask));
    end

    % Split the bounding box into two halves using a vertical line, and store each half as different masks

    bboxmask{1}(uppervert:lowervert,leftvert:leftvert+partwidth) = 1;
    bboxmask{2}(uppervert:lowervert,leftvert+partwidth+1:rightvert) = 1;

    % Determine the intersection of the original mask and each half of the bounding box to split the original mask into two halves

    for i = 1:2
        halfmask{i} = mask & bboxmask{i};
    end

    %% Load the position of islet and immune cells

    Isletcellsold = table2array(readtable([filepath,['islet ', num2str(isletno), ' glucagon position.xlsx']])); %Replace with either alpha or beta cells
    Immunecellsold = table2array(readtable([filepath,['islet ', num2str(isletno), ' cd3 position.xlsx']])); %Replace with T-, macrophage or myeloid cells

    %% Count the number of islet and immune cells inside each half-islet

    for j=1:2
    
        % Create a shape to represent the vertices of the mask representing one half of an islet

        Selection = bwboundaries(halfmask{j}); %Determine the vertices of the mask
        yMask = Selection{1}(:,1)*scaling; % Convert the x coordinates of the vertices into um
        xMask = Selection{1}(:,2)*scaling; % Convert the y coordinates of the vertices into um 
        Selection = polyshape(xMask,yMask); % Create a shape to represent the mask

        % Count the number of islet cells inside the mask

        isletcount = 0; %Initialize

        for i = 1:size(Isletcellsold,1)
            if(inpolygon(Isletcellsold(i,1),Isletcellsold(i,2),Selection.Vertices(:,1),Selection.Vertices(:,2))) %Check if the cell is inside the mask
                isletcount = isletcount + 1;
            end
        end

        isletcellcount(sno,j) = isletcount; %Store the number of islet cells in each islet half, each row corresponds to an islet and each column corresponds to one half

        % Count the number of immune cells inside the mask

        immunecount = 0; %Initialize

        for i = 1:size(Immunecellsold,1)
            if(inpolygon(Immunecellsold(i,1),Immunecellsold(i,2),Selection.Vertices(:,1),Selection.Vertices(:,2)))
                immunecount = immunecount + 1;
            end
        end

        immunecellcount(sno,j) = immunecount; %Store the number of immune cells in each islet half, each row corresponds to an islet and each column corresponds to one half
    end

end

%% Calculate the percentage of islet and immune cells in each islet half

for i = 1:size(isletcellcount,1)
    for j = 1:2
        isletcellpercent(i,j) = isletcellcount(i,j)/(isletcellcount(i,1)+isletcellcount(i,2)); %Compute the percentage of islet cells in each islet half
        immunecellpercent(i,j) = immunecellcount(i,j)/(immunecellcount(i,1)+immunecellcount(i,2)); %Compute the percentage of immune cells in each islet half
    end
end

%% Store all the percentages in a single array

isletcellall = []; %Initialize 
immunecellall = []; %Initialize

for i = 1:size(isletcellcount,1)
    for j = 1:2
        isletcellall(end+1,1) = isletcellpercent(i,j);
        immunecellall(end+1,1) = immunecellpercent(i,j);
    end
end

%% Group immune cell percentages based on islet cell percentages (bin every 5%)

binned = cell(1,20); %Create 20 bin with the first bin representing islet cell percentage of 0-5%, second cell representing 5-10% and so on until 95-100%

for i = 1:20
    binned{i} = []; %Initialize
end

for i = 1:size(isletcellall,1)
    r = floor(isletcellall(i)*100/5); %Determine the index corrsponding to a given percentage of islet cells
    r = r+1; %Factor added as the floor operator rounds down the calculations
    if r==21
        r = 20;
    end
    binned{r}(end+1,1) = immunecellall(i); %Store the percentage of immune cells in its corresponding index, based on the percentage of islet cells
end


%% Determine the mean and SEM of each bin

for i = 1:20
    meanimmune(i) = mean(binned{i}); 
    semimmune(i) = std(binned{i})/size(binned{i},1);
end

output = [meanimmune' semimmune']; %Group for easy export

%% Define a function to check if an array contains non zero values. This function returns 1 if the array contains a non zero values, otherwise it returns zero

function flag = checkforones(array) 
    flag = 0; %Initialize output of the function
    for i = 1:max(size(array,1),size(array,2)) %Iterate through the different elements of an array 
        if array(i) == 1 %Check if the element is a non zero value, and update flag if it is  
            flag = 1;
            break;
        end
    end
end
