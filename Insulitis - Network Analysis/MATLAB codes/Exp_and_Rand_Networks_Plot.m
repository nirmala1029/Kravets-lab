%% This code plots the experimental and randomized networks of an islet
%% Credit: Nirmala Balasenthilkumaran, Vira Kravets
%% Last modified: April 2024

%% Clear command window, and all saved variables

clc
clear all
close all

%% Initialization

filepath = 'C:\Nirmala\UCSD\Kravets Lab\T1D Project\All islets\'; %Initialize filepath to load files containing islet masks and the positions of various cell types in the islet

isletno = 26; %Set islet no
threshold = 20; %Distance threshold for network creation
scaling = 0.4964671; %Relationship between 1 pixel and 1 um
seed = 1; %Set the seed of the random generator

%Define the colors of cells
Fixedcellscolor = '#5D3FD3'; %Set color for plotting the fixed alpha cells
Randomcellscolor = '#D21404'; %Set color for plotting the randomized T-cells

%% Load the positions of fixed cells and random cells, and the islet masks

Randomizecellsold = table2array(readtable([filepath,['islet ', num2str(isletno), ' cd3 position.xlsx']])); %Load the positions of the cells to be randomized - T-cells
Fixedcellsold = table2array(readtable([filepath,['islet ', num2str(isletno), ' glucagon position.xlsx']])); %Load the positions of the cells to be fixed - alpha cells

maskOuter = logical(imread([filepath,['islet ', num2str(isletno), ' mask 20.tif']])); %Load the mask encompassing islets cells and immune cells, 20 um from the islet rim
isletmask = logical(imread([filepath,['islet ', num2str(isletno), ' mask 0.tif']])); %Load the mask representing the islet

immask = maskOuter - isletmask; %Create a mask to represent the periphery of an islet (up to 20um from the rim of the islet)

%% Create boundaries to represent each mask

% Create a shape to represent the vertices of the mask representing the islet and its periphery

Selection = bwboundaries(maskOuter); %Determine the vertices of the mask
yMask = Selection{1}(:,1)*scaling; %Convert the y coordinates of the vertices into um
xMask = Selection{1}(:,2)*scaling; %Convert the x coordinates of the vertices into um
Selection = polyshape(xMask,yMask); %Create a shape to represent the mask

% Create a shape to represent the vertices of the mask representing the islet 

isletselection = bwboundaries(isletmask); %Determine the vertices of the mask
yMask = isletselection{1}(:,1)*scaling; %Convert the y coordinates of the vertices into um
xMask = isletselection{1}(:,2)*scaling; %Convert the x coordinates of the vertices into um
isletselection = polyshape(xMask,yMask); %Create a shape to represent the mask

%% Generate a list of all available positions that the cells to be randomized could be situated in

% Create set of all available coordinates in the periphery of the islet

coordsOutside = []; %Initialize
[coordsOutside(:,2) coordsOutside(:,1)] = find(immask == 1); %Determine all regions where cells could be present in the mask, outside the islet
coordsOutside = coordsOutside*scaling; %Determine the coordinates of the positions

%Create set of all available coordinates inside the islet

coordsInside = []; %Initialize
[coordsInside(:,2) coordsInside(:,1)] = find(isletmask == 1); %Determine all regions where cells could be present in the mask, inside the islet
coordsInside = coordsInside*scaling; %Determine the coordinates of the positions


% Calculate distance of all available coordinates to islet rim
    
distCoordOutsidetoIslet = round(min(pdist2(coordsOutside,isletselection.Vertices),[],2)); 
distCoordInsidetoIslet = round(min(pdist2(coordsInside,isletselection.Vertices),[],2));

%% Determine the list of cells inside and outside the islet

%Initialization 

Randomcells = []; %List of cells to be randomized within 20 um of islet rim (inside mask)
con = []; %Set con = 1 if the cells to be randomized is inside islet, else set it to 0

%Determine cells within 20 um of islet rim and determine con by checking if the cell is inside the islet

for i = 1:size(Randomizecellsold,1)
    if(inpolygon(Randomizecellsold(i,1),Randomizecellsold(i,2),Selection.Vertices(:,1),Selection.Vertices(:,2))) %Check if cells are inside the mask
        Randomcells(end+1,1) = Randomizecellsold(i,1);
        Randomcells(end,2) = Randomizecellsold(i,2);
        con(end+1,1) = 0;
        if (inpolygon(Randomizecellsold(i,1),Randomizecellsold(i,2),isletselection.Vertices(:,1),isletselection.Vertices(:,2))) %Check if cells are inside the islet
            con(end,1) = 1;
        end
    end
end

% Retain only the fixed cells that are within 20 um of islet rim

Fixedcells = []; %Initialize

for i = 1:size(Fixedcellsold,1)
    if (inpolygon(Fixedcellsold(i,1),Fixedcellsold(i,2),Selection.Vertices(:,1),Selection.Vertices(:,2)))
        Fixedcells(end+1,:) = Fixedcellsold(i,:);
    end
end

% Separate the cells that are inside and outside the islet by checking con

RadomcellsOutside = [];
RandomcellsInside = [];

for i = 1:size(con)
    if con(i) == 1
        RandomcellsInside(end+1,:) = Randomcells(i,:);
    else
        RadomcellsOutside(end+1,:) = Randomcells(i,:);
    end
end

%% Construct experimental network

% Determine number of cells
nFixedcells = size(Fixedcells,1); %Compute the number of fixed cells in the islet
nRandomcells = size(Randomcells,1); %Compute the number of cells to be randomized in the islet

% Store list of cells with type in array
twocells = [Fixedcells, ones(nFixedcells,1); Randomcells, 2*ones(nRandomcells,1)];

% Define ranges for each cell type
Fixedcellsrange = 1:nFixedcells;
Randomcellsrange = nFixedcells+1:nFixedcells+nRandomcells;

ExpNetwork = zeros(size(twocells,1),size(twocells,1)); %Initialize network

%Create proximity-based network
for i =1:size(twocells,1)
    for j = 1:size(twocells,1)
        dist = distancebetweenpoints(twocells(i,1:2),twocells(j,1:2));
        if dist<threshold && twocells(i,3) ~= twocells(j,3)
            ExpNetwork(i,j) = 1;
        end
    end
end

ExpNetwork = ExpNetwork - diag(diag(ExpNetwork)); %Remove self-links
ExpNetworkgraph = graph(ExpNetwork); %Create a graph of the network

%% Plot the experimental network

figure;
ExpnetworkPlot = plot(ExpNetworkgraph,'Xdata',twocells(:,1),'YData',twocells(:,2),'NodeLabel',''); %Plot the network
highlight(ExpnetworkPlot,Fixedcellsrange,'NodeColor',Fixedcellscolor,'MarkerSize',7); %Highlight the fixed cells in a user-defined color
highlight(ExpnetworkPlot,Randomcellsrange,'NodeColor',Randomcellscolor,'MarkerSize',7); %Highlight the random-cells in a user defined color
highlight(ExpnetworkPlot,ExpNetworkgraph,'EdgeColor','#000000','LineWidth',2); %Highlight the links
xlim([50 500]); %Set x-axis scale
ylim([50 450]); %Set y-axis scale
xticks(50:100:450); %Set xticks
yticks(50:100:450); %Set yticks
ax = gca;
ax.XAxis.FontSize = 12; %Set font size of x-axis
ax.YAxis.FontSize = 12; %Set font size of y-axis
set(gca,'FontName','Arial', 'LineWidth', 2); %Bold graph and font

%% Randomization of cell positions

% Apply normal distribution to randomly generate a list of random distances between a cell and islet rim 

Randomcellsrand = []; %Initialize a variable to store randomly generated cells

%For cells outside the islet

if ~isempty(RadomcellsOutside)

    %Determine mean and standard deviation of distances of cells to islet

    distCellsOutsidetoIslet = min(pdist2(RadomcellsOutside,isletselection.Vertices),[],2); %Determime the distance between the cells and islet rim
    meandistOutside = mean(distCellsOutsidetoIslet); %Compute mean of distances
    stddistOutside = std(distCellsOutsidetoIslet); %Compute the standard deviation of distances

    %Create a normal distribution of distances

    rng(1); %Set random generator (fixed to 1)
    randdistOutside =normrnd(meandistOutside,stddistOutside,1000,1);
    randdistOutside = round(randdistOutside);

    %Randomly pick distances

    rng(seed); %Set seed of random generator

    randcount = size(RadomcellsOutside,1); %Determine the number of cells outside the islet
    indicesOutside = randi(size(randdistOutside,1),1,randcount); %Randomly pick a list of indices for distances from the normal distribution

    %Randomly assign coordinates for each random distance

    for i = 1:randcount
        dist = randdistOutside(indicesOutside(i)); %Determine the value of random distance for the selected index
        if dist > max(distCoordOutsidetoIslet)
            dist = max(distCoordOutsidetoIslet); %Assign maximum distance if the random distance is greater than the maximum available distance 
        end
        if dist < min(distCoordOutsidetoIslet)
            dist = min(distCoordOutsidetoIslet); %Assign minimum distance if the random distance is lesser than the minimum available distance 
        end

        indicesdist = find(distCoordOutsidetoIslet == dist); %Find the list of available coordinates with distance equal to the random distance
            
        %Randomly select a coordinate from the mask with distance equal to the random distance
        index = randi(size(indicesdist,1),1); 
        Randomcellsrand(end+1,:) = coordsOutside(indicesdist(index),:); 
    end

end

%For cells inside the islet

if ~isempty(RandomcellsInside)

    %Determine mean and standard deviation of distances of cells to islet

    distCellsInsidetoIslet = min(pdist2(RandomcellsInside,isletselection.Vertices),[],2); %Determime the distance between the cells and islet rim
    meandistInside = mean(distCellsInsidetoIslet); %Compute mean of distances
    stddistInside = std(distCellsInsidetoIslet); %Compute standard deviation of distances

    %Create a normal distribution of distances

    rng(1); %Set random generator (fixed to 1)
    randdistInside =normrnd(meandistInside,stddistInside,1000,1);
    randdistInside = round(randdistInside);

    %Randomly pick distances

    rng(seed); %Set seed of random generator
    randcount = size(RandomcellsInside,1); %Determine the number of cells inside the islet
    indicesInside = randi(size(randdistInside,1),1,randcount);  %Randomly pick a list of indices for distances from the normal distribution

    %Randomly assign coordinates for each random distance

    for i = 1:randcount
        dist = randdistInside(indicesInside(i)); %Determine the value of random distance for the selected index
        if dist > max(distCoordInsidetoIslet)
            dist = max(distCoordInsidetoIslet); %Assign maximum distance if the random distance is greater than the maximum available distance 
        end
        if dist < min(distCoordInsidetoIslet)
            dist = min(distCoordInsidetoIslet); %Assign minimum distance if the random distance is lesser than the minimum available distance 
        end

        indicesdist = find(distCoordInsidetoIslet == dist); %Find the list of available coordinates with distance equal to the random distance
        
        %Randomly select a coordinate from the mask with distance equal to the random distance
        index = randi(size(indicesdist,1),1); 
        Randomcellsrand(end+1,:) = coordsInside(indicesdist(index),:); 
    end
end

%% Create random network

% Determine number of cells
nFixedcells = size(Fixedcells,1); %Compute the number of cells to be fixed
nRandomcells = size(Randomcellsrand,1); %Compute the number of cells to be randomized

% Store list of cells with type in array
twocells = [Fixedcells, ones(nFixedcells,1); Randomcellsrand, 2*ones(nRandomcells,1)];

%Define cell ranges
Fixedcellsrange = 1:nFixedcells;
Randomcellsrange = nFixedcells+1:nFixedcells+nRandomcells;

RandNetwork = zeros(size(twocells,1),size(twocells,1)); %Initialize network

%Create proximity-based network
for i =1:size(twocells,1)
    for j = 1:size(twocells,1)
        dist = distancebetweenpoints(twocells(i,1:2),twocells(j,1:2));
        if dist<threshold && twocells(i,3) ~= twocells(j,3)
            RandNetwork(i,j) = 1;
        end
    end
end

RandNetwork = RandNetwork - diag(diag(RandNetwork)); %Remove self links
RandNetworkgraph = graph(RandNetwork);  %Create a graph of the random network

%% Plot random network

figure;
RandnetworkPlot = plot(RandNetworkgraph,'Xdata',twocells(:,1),'YData',twocells(:,2),'NodeLabel',''); %Plot the network
highlight(RandnetworkPlot,Fixedcellsrange,'NodeColor',Fixedcellscolor,'MarkerSize',7); %Highlight the fixed cells using a user-deifned color
highlight(RandnetworkPlot,Randomcellsrange,'NodeColor',Randomcellscolor,'MarkerSize',7); %Highlight the random cells using a user-defined color
highlight(RandnetworkPlot,RandNetworkgraph,'EdgeColor','#000000','LineWidth',2); %Highlight the links
xlim([50 500]); %Set scale for x-axis
ylim([50 450]); %Set scale for y-axis
xticks(50:100:450); %Set xticks
yticks(50:100:450); %Set yticks
ax = gca;
ax.XAxis.FontSize = 12; %Set font size of x-axis
ax.YAxis.FontSize = 12; %Set font size of y-axis
set(gca,'FontName','Arial', 'LineWidth', 2); %Bold graph and font

%% Define a function to determine the distance between two points

function dist = distancebetweenpoints(A,B)
    dist = sqrt((A(1)-B(1))^2 + (A(2)-B(2))^2); 
end
