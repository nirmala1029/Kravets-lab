%% This code randomly selects 10 islets from the list of islets classified as early-stage insulitis, and then plots the positions of alpha-, beta- and T-cells
%% Credit: Nirmala Balasenthilkumaran, Vira Kravets
%% Last modified - April 2024

%% Clear command window, and all saved variables

clc
clear all
close all

%% Initializion

filepath = 'C:\Nirmala\UCSD\Kravets Lab\T1D Project\All islets\'; %Initialize filepath to load files containing islet masks and the positions of various cell types in the islet

% Define colors for plotting cells
Bcellscolor = '#FFAC1C';
Tcellscolor = '#D21404';
Acellscolor = '#5D3FD3';

scaling = 0.4964671; %Relationship between 1 pixel and 1 um

%% Pick 10 islets at random

isletset = [58 131 72 27 34 78 81 93 26 69 126 101 92 25 97 66 53 52 7 67 132 64 90 133 56 59 55 48 114 115 62 111 47 80 127 70 96 71 89 134 31 28 54 23]; %List the islets classified as early-stage insulitis
rng(1); %Set seed
indices = randperm(size(isletset,2),10); %Randomly pick 10 indices 

%% Plot the cell positions of the randomly picked islets

for i = 1:size(indices,2) 

    isletno = isletset(indices(i)); %Obtain the islet ID from the index
    isletlist(i) = isletno; %Store the islet ID of the randomly picked early-stage islets

    % Load the positions of alpha-, beta- and T-cells in an islet, and its mask

    Tcellsold = table2array(readtable([filepath,['islet ', num2str(isletno), ' cd3 position.xlsx']])); %Load T-cell positions
    Acellsold = table2array(readtable([filepath,['islet ', num2str(isletno), ' glucagon position.xlsx']])); %Load alpha-cell positions
    Bcellsold = table2array(readtable([filepath,['islet ', num2str(isletno), ' insulin position.xlsx']])); %Load beta-cell positions

    % Create a shape to represent the vertices of the mask representing the islet and its periphery

    Mask = logical(imread([filepath,['islet ', num2str(isletno), ' mask 20.tif']])); %Load the mask encompassing islets cells and immune cells, 20 um from the islet rim 
    Selection = bwboundaries(Mask); %Determine the vertices of the mask
    yMask = Selection{1}(:,1)*scaling; %Convert the y coordinates of the vertices into um
    xMask = Selection{1}(:,2)*scaling; %Convert the x coordinates of the vertices into um
    Selection = polyshape(xMask,yMask); %Create a shape to represent the mask

    %% Retain only those cells that are inside the 20um mask

    Tcells = []; %Initialize T-cells
    Acells = []; %Initialize alpha-cells
    Bcells = []; %Initialize beta-cells

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

    %% Plot the positions of alpha,beta and T-cells

    figure;
    hold on;
    acellsplot = scatter(Acells(:,1),Acells(:,2),'filled','MarkerFaceColor',Acellscolor); %Plot alpha-cells
    bcellsplot = scatter(Bcells(:,1),Bcells(:,2),'filled','MarkerFaceColor',Bcellscolor); %Plot beta-cells
    tcellsplot = scatter(Tcells(:,1),Tcells(:,2),'filled','MarkerFaceColor',Tcellscolor); %Plot T-cells
    acellsplot.SizeData = 70; %Change size  
    bcellsplot.SizeData = 70; %Change size
    tcellsplot.SizeData = 70; %Change size
    xlim([50 300]); %Set x-axis
    ylim([50 240]); %Set y-axis
    ax = gca;
    ax.XAxis.FontSize = 12; %Set font size of x-axis
    ax.YAxis.FontSize = 12; %Set font size of y-axis 
    set(gca,'FontName','Arial','XTick',50:100:250,'YTick',50:100:250,'LineWidth', 2); %Bold graph and font
end