%% This code computes the insulitis degree for all the islets of a mouse, sorts the islets according to the degree and plots a spectrum of insulitis progression
%% Credit: Nirmala Balasenthilkumaran, Vira Kravets
%% Last modified: April 2024

%% Clear command window, and all saved variables

clc
clear all
close all

%% Initialization

filepath = 'C:\Nirmala\UCSD\Kravets Lab\T1D Project\All islets\'; %Initialize filepath to load files containing islet masks and the positions of various cell types in the islet
scaling = 0.4964671; %Relationship between 1 pixel and 1 um
colors = {'#FF00FF','#D2AFFF','#7B1FA2','#303F9F','#1976D2','#40E0D0','#388E3C','#FBC02D','#D9B99B','#FB8C00','#D32F2F'}; %Define to colors used to plot different mice
mouseno = 11; %Set mouse number of the mouse of interest
sno = 0; %Initialize looping variable

%% Calculate number of cells of each type in every islet

for isletno = 113:134 %List the islets correspondidng to the mouse of interest

    %% Import the files containing the X and Y coordinates of each T- and beta- cell

    Tcellsold = table2array(readtable([filepath,['islet ', num2str(isletno), ' cd3 position.xlsx']])); %Load the positions of T-cells
    Bcellsold = table2array(readtable([filepath,['islet ', num2str(isletno), ' insulin position.xlsx']])); %Load the positions of beta-cells

    %% Import the masks defining islet boundary and its periphery and create an outline

    maskouter = logical(imread([filepath,['islet ', num2str(isletno), ' mask 60.tif']])); %Load the mask encompassing islets cells and immune cells, 60 um from the islet rim
    outerSelection = bwboundaries(maskouter); %Determine the vertices of the mask
    yMask = outerSelection{1}(:,1)*scaling; %Convert the y coordinates of the vertices into um
    xMask = outerSelection{1}(:,2)*scaling; %Convert the x coordinates of the vertices into um
    outerSelection = polyshape(xMask,yMask); %Create a shape to represent the mask

    isletmask = logical(imread([filepath,['islet ', num2str(isletno), ' mask 0.tif']])); %Load the mask representing the islet

    %% Count the number of T-cells and beta-cells inside the mask
    
    tcellcount = 0; %Initialize

    for i = 1:size(Tcellsold,1)
        if(inpolygon(Tcellsold(i,1),Tcellsold(i,2),outerSelection.Vertices(:,1),outerSelection.Vertices(:,2)))
            tcellcount = tcellcount + 1;
        end
    end

    bcellcount = 0; %Initialize

    for i = 1:size(Bcellsold,1)
        if(inpolygon(Bcellsold(i,1),Bcellsold(i,2),outerSelection.Vertices(:,1),outerSelection.Vertices(:,2)))
            bcellcount = bcellcount + 1;
        end
    end

    sno = sno+1; %Update looping variable

    Op(sno,1) = isletno; %Store the islet number in the first column of the output file
    Op(sno,2) = tcellcount/bcellcount; %Calculate insulitis degree (ratio of T-cells and beta-cells)
    Op(sno,3) = bwarea(isletmask)*scaling*scaling; %Compute islet area

end

%% Sort the islets in ascending order based on the insulitis degree

SortedOp = sortrows(Op,2);

%% Isolate the parameters for plotting results

degree = SortedOp(:, 2); % Isolate insulitis degree of each islet
area = SortedOp(:,3); %Isolate the area of each islet

%% Plot insulitis progression

figure;
hold on;

% Create a bar graph with the x axis representing islet ID, and y axis representing insulitis degree and islet area

for i = 1:size(Op,1)
    degreebar = bar(i,degree(i),'FaceColor',colors{mouseno},'EdgeColor',colors{mouseno},'LineWidth',2); %Plot insulitis degree
    areabar = bar(i,-area(i)/10000,'FaceColor','#607D8B','EdgeColor','#607D8B','LineWidth',2); %Plot area  
end

% Define labels
set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2); %Bold font
ticks = get(gca,'YTick');
ax = gca;
ax.YTick = unique(sort([ax.YTick,-ax.YTick])); %Isolate y axis scales
ax.YTickLabel = arrayfun(@num2str,abs(ax.YTick),'UniformOutput',false); %Make negative y axis scale positive
ax.XTick = []; %Remove x-axis scale

