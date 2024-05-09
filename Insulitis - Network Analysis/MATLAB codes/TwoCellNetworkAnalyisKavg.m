%% This code creates proximity-based two-cell networks for multiple islets, and then creates an array containing average links (of the two cell network) for each islet 
%% Credit: Nirmala Balasenthilkumaran, Vira Kravets
%% Last modified - April 2024

%% Clear command window, and all saved variables

clc
clear all
close all

%% Initialization

filepath = 'C:\Nirmala\UCSD\Kravets Lab\T1D Project\All islets\'; %Initialize filepath to load files containing islet masks and the positions of various cell types in the islet

Kavgcell = []; %Initialize an array to store number of links of each cell in all islets
sno = 0; %Initialize looping variable
scaling = 0.4964671; %Relationship between 1 pixel and 1 um
threshold = 20; %Distance threshold for network creation

%% Network analysis of user selected islets

for isletno = 1:134 %Update this line with a list of islets to be analyzed

    sno = sno + 1; %Update looping variable

    %% Import coordinates of cells

    Mcellsold = table2array(readtable([filepath1,['islet ', num2str(isletno), ' cd11c position.xlsx']])); %Import coordinates of myeloid cells
    Tcellsold = table2array(readtable([filepath1,['islet ', num2str(isletno), ' cd3 position.xlsx']])); %Import coordinates of T-cells
    Macrophagesold = table2array(readtable([filepath1,['islet ', num2str(isletno), ' f480 position.xlsx']])); %Import coordinates of macrophages
    Acellsold = table2array(readtable([filepath1,['islet ', num2str(isletno), ' glucagon position.xlsx']])); %Import coordinates of alpha cells
    Bcellsold = table2array(readtable([filepath1,['islet ', num2str(isletno), ' insulin position.xlsx']])); %Import coordinates of beta cells

    %Load mask and create selection

    Mask = logical(imread([filepath2,['islet ', num2str(isletno), ' mask 20.tif']])); %Load the mask encompassing islets cells and immune cells, 20 um from the islet rim
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
    
    %% Construct the two-cell network
    
    %Define the cell types to be analyzed (averages wrt the first cell type)

    celltype1 = Mcells; 
    celltype2 = Macrophages;

    %Determine the number of cells of selected cell types

    ncelltype1 = size(celltype1,1);
    ncelltype2 = size(celltype2,1);

    % Store list of cells with coordinates as first 2 columns and the type as third column in an array
    twocells = [celltype1, ones(ncelltype1,1); celltype2, 2*ones(ncelltype2,1)];

    %Define ranges of cell type
    Celltype1range = 1:ncelltype1;
    Celltype2range = ncelltype1+1:ncelltype1+ncelltype2;

    TwoCellsNetwork = zeros(size(twocells,1),size(twocells,1)); %Initialize the network

    %Create a two cell network based on proximity

    for i = 1:size(twocells,1)
        for j = 1:size(twocells,1)
            dist = distancebetweenpoints(twocells(i,1:2),twocells(j,1:2));
            if dist<threshold && twocells(i,3) ~= twocells(j,3)
                TwoCellsNetwork(i,j) = 1;
            end
        end
    end

    TwoCellsNetwork = TwoCellsNetwork - diag(diag(TwoCellsNetwork)); %Remove self-links

    %% Compute Kavg

    Kavgcount = []; %Define an empty array to store the number of links of each cell in an islet

    for i = 1:size(twocells,1) 
        count = 0; %Initalize number of links of a cell to zero
        % Count the number of links of a cell
        for j = 1:size(twocells,1)
            if TwoCellsNetwork(i,j) == 1
                count = count + 1;
            end
        end
        Kavgcount(end+1) = count; %Store the number of links of a cell
        Kavgcell(end+1,1) = count; %Store the number of links of a cell 
    end

    Kavgislet(sno,1) = mean(Kavgcount); %Determine the average number of links of the network (wrt to one cell type)

end

%% Define a function to calculate the distance between two points

function dist = distancebetweenpoints(A,B)
    dist = sqrt((A(1)-B(1))^2 + (A(2)-B(2))^2); 
end
