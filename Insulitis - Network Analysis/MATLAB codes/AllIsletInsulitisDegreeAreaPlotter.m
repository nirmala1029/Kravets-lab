%% This code computes the insulitis degree and t-cell density for all available islets, sorts the islets according to the degree and plots a spectrum of insulitis progression
%% Credit: Nirmala Balasenthilkumaran, Vira Kravets
%% Last modified - April 2024

%% Clear command window, and all saved variables

clc
clear all
close all

%% Initializion

filepath = 'C:\Nirmala\UCSD\Kravets Lab\T1D Project\All islets\'; %Initialize filepath to load files containing islet masks and the positions of various cell types in the islet
scaling = 0.4964671; %Relationship between 1 pixel and 1 um
sno = 0; %Initialize looping variable

%% Count the number of T-cells and beta cells and compute insulitis degree in each islet

for isletno = 1:134 %Update this line with a list of islets to be analyzed

    %% Import the files containing the X and Y coordinates of each T- and beta- cell

    Tcellsold = table2array(readtable([filepath,['islet ', num2str(isletno), ' cd3 position.xlsx']])); %Load the positions of T-cells
    Bcellsold = table2array(readtable([filepath,['islet ', num2str(isletno), ' insulin position.xlsx']])); %Load the positions of beta cells

    %% Import the masks defining islet boundary, and its periphery

    isletmask = logical(imread([filepath,['islet ', num2str(isletno), ' mask 0.tif']])); %Load the mask representing the islet
    maskouter = logical(imread([filepath,['islet ', num2str(isletno), ' mask 60.tif']])); %Load the mask encompassing islets cells and immune cells, 60 um from the islet rim
    
    % Create a shape to represent the vertices of the mask representing the islet and its periphery

    outerSelection = bwboundaries(maskouter); %Determine the vertices of the mask
    yMask = outerSelection{1}(:,1)*scaling; %Convert the y coordinates of the vertices into um
    xMask = outerSelection{1}(:,2)*scaling; %Convert the x coordinates of the vertices into um
    outerSelection = polyshape(xMask,yMask); %Create a shape to represent the mask

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

    maskperiph = maskouter - isletmask; %Compute the regions of the periphery of the islet

    Op(sno,1) = isletno; %Store the islet number in the first column of the output file
    Op(sno,2) = tcellcount/bcellcount; %Calculate insulitis degree (ratio of T-cells and beta-cells)
    Op(sno,3) = log(Op(sno,2)); %Calculate log of insulitis degree
    Op(sno,4) = bwarea(isletmask)*scaling*scaling; %Calculate islet area
    Op(sno,5) = tcellcount; %Compute the no of T-cells
    Op(sno,6) = bcellcount; %Compute the no of beta-cells
    Op(sno,7) = tcellcount/(bwarea(maskperiph)*scaling*scaling*10^-6); %Compute T-cell density

end

%% Set Mouse number as the fifth column 

for i = 1:2
    Op(i,5) = 1;
end

for i = 3:13
    Op(i,5) = 2;
end

for i = 14:18
    Op(i,5) = 3;
end

for i = 19:24
    Op(i,5) = 4;
end

for i = 25:38
    Op(i,5) = 5;
end

for i = 39:51
    Op(i,5) = 6;
end

for i = 52:79
    Op(i,5) = 7;
end

for i = 80:98
    Op(i,5) = 8;
end

for i = 99:109
    Op(i,5) = 9;
end

for i = 110:112
    Op(i,5) = 10;
end

for i = 113:134
    Op(i,5) = 11;
end

%% Sort the islets in ascending order based on the insulitis degree

SortedOp = sortrows(Op,2);

%% Isolate the parameters for plotting results

degree = SortedOp(:, 2); %Isolate insulitis degree 
area = SortedOp(:,4); %Isolate islet area
mouseNumber = SortedOp(:,5); %Isolate mouse number 

%% Initialize legends and the colors used to represent islets of different mice

MiceNos = {'Mouse 1','Mouse 2','Mouse 3','Mouse 4','Mouse 5','Mouse 6','Mouse 7','Mouse 8', 'Mouse 9', 'Mouse 10', 'Mouse 11'}; %Set names of mouse numbers for legend
colors = {'#FF00FF','#D2AFFF','#7B1FA2','#303F9F','#1976D2','#40E0D0','#388E3C','#FBC02D','#D9B99B','#FB8C00','#D32F2F'}; %Set colors for different mice

%% Plot progression of T1D

figure;
hold on;

%Create a bar graph with the x axis representing islet ID and y axis representing insulitis degree

for i = 1:size(Op,1)
    insulitisdegreebars = bar(i,degree(i),'FaceColor',colors{mouseNumber(i)},'EdgeColor',colors{mouseNumber(i)},'LineWidth',2); %Plot insulitis degree
    areabars = bar(i,-area(i)/10000,'FaceColor','#808080','EdgeColor','#808080','LineWidth',2); %Plot area
end

%Create an empty plot to obtain a legend to represent the different colors used for each mouse

for i=1:11
    t(i) = bar(NaN,NaN,'FaceColor',colors{i},'EdgeColor',colors{i},'LineWidth',2);
end

% Define labels
set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2); %Bold font
legend(t,MiceNos,'FontSize',14,'NumColumns',5); %Define legend to highlight different mice
ax = gca;
ax.YTick = unique(sort([ax.YTick,-ax.YTick])); %Isolate y axis scales
ax.YTickLabel = arrayfun(@num2str,abs(ax.YTick),'UniformOutput',false); %Make negative y axis scale positive
ax.XTick = []; %Remove x axis scale

