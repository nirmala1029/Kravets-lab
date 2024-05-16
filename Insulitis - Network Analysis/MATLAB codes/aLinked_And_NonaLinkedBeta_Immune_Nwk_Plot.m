%% This code identifies alpha-linked and non-alpha-linked beta cells and creates proximity-based networks of these cells with immune cells
%% Credit: Nirmala Balasenthilkumaran, Vira Kravets
%% Last modified: April 2024

%% Clear command window, and all saved variables

clc
clear all
close all

%% Initialization

filepath = 'C:\Nirmala\UCSD\Kravets Lab\T1D Project\All islets\'; %Initialize filepath to load files containing islet masks and the positions of various cell types in the islet

% Define the colors of cells

Bcellscolor = '#FFAC1C'; %Set color for plotting beta-cells
Immunecellscolor = '#D21404'; %Set color for plotting immune cells
Acellscolor = '#5D3FD3'; %Set color for plotting alpha-cells

isletno = 27; %Set islet no
threshold = 20; %Distance threshold for network creation
scaling = 0.4964671; %Relationship between 1 pixel and 1 um

%% Create boundaries to represent each mask

% Create a shape to represent the vertices of the mask representing the islet 

isletmask = logical(imread([filepath,['islet ', num2str(isletno), ' mask 0.tif']])); %Load the mask representing the islet
isletSelection = bwboundaries(isletmask); %Determine the vertices of the mask
yMask = isletSelection{1}(:,1)*scaling; %Convert the y coordinates of the vertices into um
xMask = isletSelection{1}(:,2)*scaling; %Convert the x coordinates of the vertices into um
isletSelection = polyshape(xMask,yMask); %Create a shape to represent the mask

% Create a shape to represent the vertices of the mask representing the inner core of the islet (20 um inside the islet rim)

innermask = logical(imread([filepath,['islet ', num2str(isletno), ' mask -20.tif']])); %Load the mask representing the inner core of the islet
innerSelection = bwboundaries(innermask); %Determine the vertices of the mask
yMask = innerSelection{1}(:,1)*scaling; %Convert the y coordinates of the vertices into um
xMask = innerSelection{1}(:,2)*scaling; %Convert the x coordinates of the vertices into um
innerSelection = polyshape(xMask,yMask); %Create a shape to represent the mask

% Create a shape to represent the vertices of the mask representing the islet and its periphery

outermask = logical(imread([filepath,['islet ', num2str(isletno), ' mask 20.tif']])); %Load the mask encompassing islets cells and immune cells, 20 um from the islet rim
outerSelection = bwboundaries(outermask); %Determine the vertices of the mask
yMask = outerSelection{1}(:,1)*scaling; %Convert the y coordinates of the vertices into um
xMask = outerSelection{1}(:,2)*scaling; %Convert the x coordinates of the vertices into um
outerSelection = polyshape(xMask,yMask); %Create a shape to represent the mask

%% Load the positions of alpha-cells, beta-cells and immune-cells

Acellsold = table2array(readtable([filepath,['islet ', num2str(isletno), ' glucagon position.xlsx']])); %Load the positions of alpha-cells
Bcellsold = table2array(readtable([filepath,['islet ', num2str(isletno), ' insulin position.xlsx']])); %Load the positions of beta-cells
Immunecellsold = table2array(readtable([filepath,['islet ', num2str(isletno), ' cd3 position.xlsx']])); %Load the positions of immune-cells 
    
% Retain only those alpha and beta cells that are in the outer boundary of the islet (not in the inner core)  

Bcells = []; %Initialize

for i = 1:size(Bcellsold,1)
    if(inpolygon(Bcellsold(i,1),Bcellsold(i,2),outerSelection.Vertices(:,1),outerSelection.Vertices(:,2))) %Check if cells are inside the 20 um mask
        if(~inpolygon(Bcellsold(i,1),Bcellsold(i,2),innerSelection.Vertices(:,1),innerSelection.Vertices(:,2))) %Check whether cells are not a part of the inner core of the islet
            Bcells(end+1,1) = Bcellsold(i,1);
            Bcells(end,2) = Bcellsold(i,2);
        end
    end
end

Acells = []; %Initialize

for i = 1:size(Acellsold,1)
    if(inpolygon(Acellsold(i,1),Acellsold(i,2),outerSelection.Vertices(:,1),outerSelection.Vertices(:,2))) %Check if cells are inside the 20 um mask
        if(~inpolygon(Acellsold(i,1),Acellsold(i,2),innerSelection.Vertices(:,1),innerSelection.Vertices(:,2))) %Check whether cells are not a part of the inner core of the islet
            Acells(end+1,1) = Acellsold(i,1);
            Acells(end,2) = Acellsold(i,2);
        end
    end
end

% Retain only those immune cells that are within 20 um of islet

Immunecells = [];

for i = 1:size(Immunecellsold,1)
    if(inpolygon(Immunecellsold(i,1),Immunecellsold(i,2),outerSelection.Vertices(:,1),outerSelection.Vertices(:,2))) %Check if cells within the 20 um mask
        Immunecells(end+1,1) = Immunecellsold(i,1);
        Immunecells(end,2) = Immunecellsold(i,2);
    end
end


%% Create a three-cell proximity based network to identify alpha-linked beta cells

% Determine the number of cells

nImmunecells = size(Immunecells,1); %Determine the number of immune cells
nAcells = size(Acells,1); %Determine the number of alpha-cells
nBcells = size(Bcells,1); %Determine the number of beta-cells

% Store list of cells with type in array
threecells = [Immunecells, ones(nImmunecells,1); Acells, 2*ones(nAcells,1); Bcells, 3*ones(nBcells,1)];

% Define cell ranges

ImmunecellsRange = 1:nImmunecells;
AcellsRange = nImmunecells+1:nImmunecells+nAcells;    
BcellsRange = nImmunecells+nAcells+1:nImmunecells+nAcells+nBcells;

ThreeCellsNetwork = zeros(size(threecells,1),size(threecells,1)); %Initialize network

% Create the network based on proximity

for i = 1:size(threecells,1)
    for j = 1:size(threecells,1)
        dist = distancebetweenpoints(threecells(i,1:2),threecells(j,1:2));
        if dist < threshold
            ThreeCellsNetwork(i,j) = 1;
        end
    end
end

ThreeCellsNetwork = ThreeCellsNetwork - diag(diag(ThreeCellsNetwork)); %Remove self-links

%% Identify alpha-linked beta cells by checking if a beta cell is linked to an alpha cell

AlphaLinkedBeta = []; %Initialize
for i = BcellsRange
    for j = AcellsRange
        if ThreeCellsNetwork(i,j) == 1 %Check if beta cell is linked to alpha-cells
            AlphaLinkedBeta(end+1) = i; %Store the indices of the alpha-linked beta-cells
            break;
        end
    end
end

%% Identify non-alpha linked beta cells by checking if a beta cell is linked to an alpha cell

NonAlphaLinkedBeta = []; %Initialize
for i = BcellsRange
    con = 0;
    for j = AcellsRange
        if ThreeCellsNetwork(i,j) == 1 %Check if beta cell is linked to alpha-cells
            con = 1; %Set a control variable to be 1 if the beta cell is alpha-linked
        end
    end
    if con == 0 
        NonAlphaLinkedBeta(end+1) = i; %Store the index of the beta cell if it is not linked to an alpha-cells
    end
end

%% Plot the alpha and beta cells in an islet

figure;
hold on;

AcellsAllPlot = scatter(Acellsold(:,1),Acellsold(:,2),'filled','MarkerFaceColor',Acellscolor); %Plot all alpha-cells
BcellsAllPlot = scatter(Bcellsold(:,1),Bcellsold(:,2),'filled','MarkerFaceColor',Bcellscolor); %Plot all beta-cells 

AcellsAllPlot.SizeData = 70; %Change size
BcellsAllPlot.SizeData = 70; %Change size

xlim([90 300]); %Set x-axis scale
ylim([90 250]); %Set y-axis scale

ax = gca;
ax.XAxis.FontSize = 12; %Set font size of x-axis
ax.YAxis.FontSize = 12; %Set font size of y-axis
set(gca,'FontName','Arial', 'LineWidth', 2,'XTick',50:100:300,'YTick',50:100:250); %Bold graph and font

%% Plot all islet cells and highlight the alpha and beta cells in the outer ring

figure;
hold on;

AcellAllPlot = scatter(Acellsold(:,1),Acellsold(:,2),'MarkerEdgeColor',Acellscolor); %Plot all alpha-cells
BcellAllPlot = scatter(Bcellsold(:,1),Bcellsold(:,2),'MarkerEdgeColor',Bcellscolor); %Plot all beta-cells
AcellOuterPlot = scatter(Acells(:,1),Acells(:,2),'filled','MarkerFaceColor',Acellscolor); %Plot outer alpha-cells
BcellOuterPlot = scatter(Bcells(:,1),Bcells(:,2),'filled','MarkerFaceColor',Bcellscolor); %Plot outer-beta cells
    
AcellAllPlot.SizeData = 70; %Change size
BcellAllPlot.SizeData = 70; %Change size
AcellOuterPlot.SizeData = 70; %Change size
BcellOuterPlot.SizeData = 70; %Change size

xlim([90 300]); %Set x-axis scale
ylim([90 250]); %Set y-axis scale

ax = gca;
ax.XAxis.FontSize = 12; %Set font size of x-axis
ax.YAxis.FontSize = 12; %Set font size of y-axis
set(gca,'FontName','Arial', 'LineWidth', 2,'XTick',50:100:300,'YTick',50:100:250); %Bold graph and font and set xticks and yticks

%% Differentiate alpha linked beta cells and non alpha linked beta cells in the outer ring of the islet

figure;
hold on;
    
for i = AlphaLinkedBeta
    AlinkedBcellPlot = scatter(threecells(i,1),threecells(i,2),'filled','MarkerFaceColor','#FF00FF'); %Plot outer alpha-linked beta-cells
    AlinkedBcellPlot.SizeData = 70; %Change size
end
for i = NonAlphaLinkedBeta
    NAlinkedBcellPlot = scatter(threecells(i,1),threecells(i,2),'filled','MarkerFaceColor',Bcellscolor); %Plot outer non-alpha-linked beta-cells
    NAlinkedBcellPlot.SizeData = 70; %Change size
end
for i = AcellsRange
    AcellPlot = scatter(threecells(i,1),threecells(i,2),'filled','MarkerFaceColor',Acellscolor); %Plot outer-alpha cells
    AcellPlot.SizeData = 70; %Change size
end

xlim([90 300]); %Set x-axis scale
ylim([90 250]); %Set y-axis scale

ax = gca;
ax.XAxis.FontSize = 12; %Set font size of x-axis
ax.YAxis.FontSize = 12; %Set font size of y-axis
set(gca,'FontName','Arial', 'LineWidth', 2,'XTick',50:100:300,'YTick',50:100:250); %Bold graph and font

%% Construct two-cell network between alpha-linked beta cells and immune cells
    
twocells = []; %Initialize array to store cell positions

for i = AlphaLinkedBeta
    twocells(end+1,:) = threecells(i,1:2); %Copy the coordinates of alpha linked beta cells
end

nALBcells = size(twocells,1); %Compute the number of alpha linked beta cells
nImmunecells = size(Immunecells,1); %Compute number of immune cells

twocells = [twocells ones(nALBcells,1); Immunecells 2*ones(nImmunecells,1)]; % Store list of cells with type in array

%Define ranges for each cell type
AlphaLinkedBeta = 1:nALBcells;
ImmunecellsRange = nALBcells+1:nALBcells+nImmunecells;

TwoCellsNetwork1 = zeros(size(twocells,1),size(twocells,1)); %Initialize the network

%Create a two cell network based on proximity
for i =1:size(twocells,1)
    for j = 1:size(twocells,1)
        dist = distancebetweenpoints(twocells(i,1:2),twocells(j,1:2));
        if dist<threshold && twocells(i,3) ~= twocells(j,3)
            TwoCellsNetwork1(i,j) = 1;
        end
    end
end
TwoCellsNetwork1 = TwoCellsNetwork1 - diag(diag(TwoCellsNetwork1)); %Remove self-links
TwoCellsNetwork1graph = graph(TwoCellsNetwork1); %Create a graph of the network

%% Plot alpha-linked beta cell - immune cell network

figure;
ALBGraphPlot = plot(TwoCellsNetwork1graph,'Xdata',twocells(:,1),'YData',twocells(:,2),'NodeLabel',''); %Plot the network
highlight(ALBGraphPlot,AlphaLinkedBeta,'NodeColor','#FF00FF','MarkerSize',7); %Highlight the alpha-linked beta-cells in magenta
highlight(ALBGraphPlot,ImmunecellsRange,'NodeColor',Immunecellscolor,'MarkerSize',7); %Highlight the immune cells in the user defined color
highlight(ALBGraphPlot,TwoCellsNetwork1graph,'EdgeColor','#000000','LineWidth',2); %Highlight the links
xlim([90 300]); %Set x-axis scale
ylim([90 250]); %Set y-axis scale
ax = gca;
ax.XAxis.FontSize = 12; %Set font size of x-axis 
ax.YAxis.FontSize = 12; %Set font size of y-axis
set(gca,'FontName','Arial', 'LineWidth', 2,'XTick',50:100:250,'YTick',50:100:250); %Bold graph and font

%% Construct two-cell network between non alpha-linked beta cells and immune cells 

twocells = []; %Initialize array to store cell positions

for i = NonAlphaLinkedBeta
    twocells(end+1,:) = threecells(i,1:2); %Copy the coordinates of non-alphalinekd beta cells
end

nNALBcells = size(twocells,1); %Compute number of non-alpha linked beta cells
nImmunecells = size(Immunecells,1); %Compute number of immune cells

twocells = [twocells ones(nNALBcells,1); Immunecells 2*ones(nImmunecells,1)]; %% Store list of cells with type in array

%Define ranges for each cell type
NonAlphaLinkedBeta = 1:nNALBcells;
ImmunecellsRange = nNALBcells+1:nNALBcells+nImmunecells;

TwoCellsNetwork2 = zeros(size(twocells,1),size(twocells,1)); %Initialize network

%Create a proximity-based network
for i =1:size(twocells,1)
    for j = 1:size(twocells,1)
        dist = distancebetweenpoints(twocells(i,1:2),twocells(j,1:2));
        if dist<threshold && twocells(i,3) ~= twocells(j,3)
            TwoCellsNetwork2(i,j) = 1;
        end
    end
end

TwoCellsNetwork2 = TwoCellsNetwork2 - diag(diag(TwoCellsNetwork2)); %Remove self-links
TwoCellsNetwork2graph = graph(TwoCellsNetwork2); %Create a graph of the network

%% Plot non alpha-linked beta-cell - immune-cell network

figure;
NALBGraphPlot = plot(TwoCellsNetwork2graph,'Xdata',twocells(:,1),'YData',twocells(:,2),'NodeLabel',''); %Plot the network
highlight(NALBGraphPlot,NonAlphaLinkedBeta,'NodeColor',Bcellscolor,'MarkerSize',7); %Highlight the non-alpha linked beta-cells in user defined color
highlight(NALBGraphPlot,ImmunecellsRange,'NodeColor',Immunecellscolor,'MarkerSize',7); %Highlight the immune cells in user defined color
highlight(NALBGraphPlot,TwoCellsNetwork2graph,'EdgeColor','#000000','LineWidth',2); %Highlight the links
xlim([90 300]); %Set x-axis scale
ylim([90 250]); %Set y-axis scale
ax = gca;
ax.XAxis.FontSize = 12; %Set font size of x-axis
ax.YAxis.FontSize = 12; %Set font size of y-axis
set(gca,'FontName','Arial', 'LineWidth', 2,'XTick',50:100:250,'YTick',50:100:250); %Bold graph and font and set xticks and yticks


%% Define a function to determine the distance between two points

function dist = distancebetweenpoints(A,B)
    dist = sqrt((A(1)-B(1))^2 + (A(2)-B(2))^2); 
end
