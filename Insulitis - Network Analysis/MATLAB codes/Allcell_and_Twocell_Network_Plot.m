%% This code creates a network using all cell types, and also creates a network using any two user selected cell types
%% Credit: Nirmala Balasenthilkumaran, Vira Kravets
%% Last modified - April 2024

%% Clear command window, and all saved variables

clc
clear all 
close all

%% Initialization  

filepath = 'C:\Nirmala\UCSD\Kravets Lab\T1D Project\All islets\'; %Initialize filepath to load files containing islet masks and the positions of various cell types in the islet

isletno = 5; %Set islet number of choice
threshold = 20; %Set threshold for network creation
scaling = 0.4964671; %Relationship between 1 pixel and 1 um

Image = imread(['islet ', num2str(isletno), ' cd3.tif']); %Read an image to obtain size for white background

%% Load cell positions of all cell types

Mcellsold = table2array(readtable([filepath,['islet ', num2str(isletno), ' cd11c position.xlsx']])); %Load positions of myeloid cells
Tcellsold = table2array(readtable([filepath,['islet ', num2str(isletno), ' cd3 position.xlsx']])); %Load positions of T-cells
Macrophagesold = table2array(readtable([filepath,['islet ', num2str(isletno), ' f480 position.xlsx']])); %Load positions of macrophages
Acellsold = table2array(readtable([filepath,['islet ', num2str(isletno), ' glucagon position.xlsx']])); %Load positions of alpha-cells
Bcellsold = table2array(readtable([filepath,['islet ', num2str(isletno), ' insulin position.xlsx']])); %Load positions of beta-cells

%% Load mask and create a shape for selection

Mask = logical(imread([filepath,['islet ', num2str(isletno), ' mask 20.tif']])); %Load the mask encompassing islets cells and immune cells, 20 um from the islet rim
Selection = bwboundaries(Mask);
yMask = Selection{1}(:,1)*scaling; % Convert the y coordinates of the vertices into um
xMask = Selection{1}(:,2)*scaling; % Convert the x coordinates of the vertices into um
Selection = polyshape(xMask,yMask); % Create a shape to represent the mask

%% Retain only those cells that are inside the 20um mask 

% Initialize empty arrays for each cell type

Mcells = []; %Myeloid cells
Tcells = []; %T-cells
Macrophages = []; %Macrophages
Acells = []; %Alpha cells
Bcells = []; %Beta cells

for i = 1:size(Tcellsold,1)
    if(inpolygon(Tcellsold(i,1),Tcellsold(i,2),Selection.Vertices(:,1),Selection.Vertices(:,2)));
        Tcells(end+1,1) = Tcellsold(i,1);
        Tcells(end,2) = Tcellsold(i,2);
    end
end

for i = 1:size(Acellsold,1)
    if(inpolygon(Acellsold(i,1),Acellsold(i,2),Selection.Vertices(:,1),Selection.Vertices(:,2)))
        Acells(end+1,1) = Acellsold(i,1);
        Acells(end,2) = Acellsold(i,2);
    end
end

for i = 1:size(Bcellsold,1)
    if(inpolygon(Bcellsold(i,1),Bcellsold(i,2),Selection.Vertices(:,1),Selection.Vertices(:,2)));
        Bcells(end+1,1) = Bcellsold(i,1);
        Bcells(end,2) = Bcellsold(i,2);
    end
end    

for i = 1:size(Mcellsold,1)
    if(inpolygon(Mcellsold(i,1),Mcellsold(i,2),Selection.Vertices(:,1),Selection.Vertices(:,2)));
        Mcells(end+1,1) = Mcellsold(i,1);
        Mcells(end,2) = Mcellsold(i,2);
    end
end  

for i = 1:size(Macrophagesold,1)
    if(inpolygon(Macrophagesold(i,1),Macrophagesold(i,2),Selection.Vertices(:,1),Selection.Vertices(:,2)));
        Macrophages(end+1,1) = Macrophagesold(i,1);
        Macrophages(end,2) = Macrophagesold(i,2);
    end
end

Imagep = ones(size(Image,1),size(Image,2)); %Create a white background for the network image

%% Construct all-cell network

% Determine number of cells of each cell type

nMcells = size(Mcells,1);
nTcells = size(Tcells,1);
nMacrophages = size(Macrophages,1);
nAcells = size(Acells,1);
nBcells = size(Bcells,1);

%Create a list of all cells with coordinates in first 2 columns and type as third column

allcells = [Mcells, ones(nMcells,1); Tcells, 2*ones(nTcells,1); Macrophages, 3*ones(nMacrophages,1); Acells, 4*ones(nAcells,1); Bcells, 5*ones(nBcells,1)];

%Define cell ranges

McellsRange = 1:nMcells;
TcellsRange = nMcells+1:nMcells+nTcells;
MacrophagesRange = nMcells+nTcells+1:nMcells+nTcells+nMacrophages;
AcellsRange = nMcells+nTcells+nMacrophages+1:nMcells+nTcells+nMacrophages+nAcells;
BcellsRange = nMcells+nTcells+nMacrophages+nAcells+1:nMcells+nTcells+nMacrophages+nAcells+nBcells;

AllCellsNetwork = zeros(size(allcells,1),size(allcells,1)); %Initialize network

%Create all cell network based on proximity

for i =1:size(allcells,1)
    for j = 1:size(allcells,1)
        dist = distancebetweenpoints(allcells(i,1:2),allcells(j,1:2));
        if dist<threshold 
            AllCellsNetwork(i,j) = 1;
        end
    end
end

AllCellsNetwork = AllCellsNetwork - diag(diag(AllCellsNetwork)); %Remove self linking
AllCellsNetworkgraph = graph(AllCellsNetwork); %Create a graph of the network

%% Plot all-cell network

%Initialize colors for each cell type
mcellscolor = '#32CD32';
tcellscolor = '#D21404';
macrophagescolor = '#00A3A3';
acellscolor = '#5D3FD3';
bcellscolor = '#FFAC1C';

figure;
imshow(Imagep); %Display background image
hold on;
p = plot(AllCellsNetworkgraph,'Xdata',allcells(:,1)/scaling,'YData',allcells(:,2)/scaling); %Plot the network
set(gca, 'YDir','reverse'); %Reverse axes to ensure orientation

%Increase the size of data points and change the color of each cell based on its cell type
highlight(p,McellsRange,'NodeColor',mcellscolor,'MarkerSize',7); 
highlight(p,TcellsRange,'NodeColor',tcellscolor,'MarkerSize',7);
highlight(p,MacrophagesRange,'NodeColor',macrophagescolor,'MarkerSize',7);
highlight(p,AcellsRange,'NodeColor',acellscolor,'MarkerSize',7);
highlight(p,BcellsRange,'NodeColor',bcellscolor,'MarkerSize',7);

highlight(p,AllCellsNetworkgraph,'EdgeColor','#000000','LineWidth',2); %Increase size of links
grid on;
set(gca,'FontName','Open Sans','FontSize',12,'FontWeight','Bold', 'LineWidth', 2); %Bold graph and font

%% Construct two-cell network

% Set cell types

celltype1 = Mcells;
celltype2 = Tcells;

% Determine number of cells of each type
ncelltype1 = size(celltype1,1);
ncelltype2 = size(celltype2,1);

% Store list of cells with coordinates as first 2 columns and the type as third column in an array
twocells = [celltype1, ones(ncelltype1,1); celltype2, 2*ones(ncelltype2,1)];

% Define cell ranges
Celltype1Range = 1:ncelltype1;
Celltype2Range = ncelltype1+1:ncelltype1+ncelltype2;

TwoCellsNetwork = zeros(size(twocells,1),size(twocells,1)); %Initialize network

%Create a two cell network based on proximity

for i =1:size(twocells,1)
    for j = 1:size(twocells,1)
        dist = distancebetweenpoints(twocells(i,1:2),twocells(j,1:2));
        if dist<threshold && twocells(i,3) ~= twocells(j,3)
            TwoCellsNetwork(i,j) = 1;
        end
    end
end

TwoCellsNetwork = TwoCellsNetwork - diag(diag(TwoCellsNetwork)); %Remove self-links
TwoCellsNetworkgraph = graph(TwoCellsNetwork); %Create a graph of the network

%% Plot two-cell network

celltype1color= mcellscolor;
celltype2color = tcellscolor;

figure;
imshow(Imagep); %Display background image
hold on;
p = plot(TwoCellsNetworkgraph,'Xdata',twocells(:,1)/scaling,'YData',twocells(:,2)/scaling); %Plot the network
set(gca, 'YDir','reverse'); %Reverse axes to ensure orienttion

%Increase the size of data points and change the color of each cell based on its cell type
highlight(p,Celltype1Range,'NodeColor',celltype1color,'MarkerSize',7);
highlight(p,Celltype2Range,'NodeColor',celltype2color,'MarkerSize',7);

highlight(p,TwoCellsNetworkgraph,'EdgeColor','#000000','LineWidth',2); %Increase size of links
grid on;
set(gca,'FontName','Open Sans','FontSize',12,'FontWeight','Bold', 'LineWidth', 2); %Bold graph and font

%% Define a function to calculate the distance between two points

function dist = distancebetweenpoints(A,B)
    dist = sqrt((A(1)-B(1))^2 + (A(2)-B(2))^2); 
end
