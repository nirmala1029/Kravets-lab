%% This code creates experimental and randomized proximity-based networks of an islet, and computes average links of these networks
%% Credit: Nirmala Balasenthilkumaran, Kravets lab

clc
clear all
close all

%% Initialize various parameters

filepath = 'D:\Nirmala\Insulitis project\All islets\'; %Update filepath 
threshold = 15; %Distance threshold for network creation
scaling = 0.4964674; %Relationship between 1 pixel and 1 um
numseeds = 100; %Set the number of seeds

%Initialize variables outside loop
Kavgcellexp = [];
Kavgrandcell = cell(1,numseeds);
Kavgisletexp = [];
Kavgrandislet = cell(1,numseeds);
    
for isletno = 1:134 %List the islets to be analyzed
    
    %% Load the positions of alpha cells and immune-cells, and the islet mask
    
    Tcellsold = table2array(readtable([filepath,['islet ', num2str(isletno), ' cd11c position.xlsx']])); %Can load other immune cells changing the input file of this variable
    Acellsold = table2array(readtable([filepath,['islet ', num2str(isletno), ' glucagon position.xlsx']]));

    Mask = logical(imread([filepath,['islet ', num2str(isletno), ' mask 20.tif']]));
    isletmask = logical(imread([filepath,['islet ', num2str(isletno), ' mask 0.tif']]));
    immask = Mask - isletmask; %Create a mask to represent the periphery of an islet (up to 20um from the rim of the islet)

    %Create a boundary to represent each mask

    Selection = bwboundaries(Mask);
    yMask = Selection{1}(:,1)*scaling;
    xMask = Selection{1}(:,2)*scaling;
    Selection = polyshape(xMask,yMask);

    isletselection = bwboundaries(isletmask);
    yMask = isletselection{1}(:,1)*scaling;
    xMask = isletselection{1}(:,2)*scaling;
    isletselection = polyshape(xMask,yMask);

    %% Generate a list of all available positions that the immune cells could be situated in

    %Create set of all available coordinates in the periphery of the islet

    coordsOutside = [];
    [coordsOutside(:,2) coordsOutside(:,1)] = find(immask == 1);
    coordsOutside = coordsOutside*scaling;

    %Create set of all available coordinates inside the islet

    coordsInside = [];
    [coordsInside(:,2) coordsInside(:,1)] = find(isletmask == 1);
    coordsInside = coordsInside*scaling;

    %Calculate distance of all available coordinates to islet rim
    
    distCoordOutsidetoIslet = min(pdist2(coordsOutside,isletselection.Vertices),[],2); 
    distCoordInsidetoIslet = min(pdist2(coordsInside,isletselection.Vertices),[],2);
    distCoordInsidetoIslet = round(distCoordInsidetoIslet);
    distCoordOutsidetoIslet = round(distCoordOutsidetoIslet);

    %% Determine the list of immune cells inside and outside the islet 

    %Initialize variables 

    Tcells = []; %List of immune cells within 20 um of islet rim (inside mask)
    con = []; %Set con = 1 if immune cell is inside islet, else set it to 0

    %Determine immune cells within 20 um of islet rim and determine con by checking if the cell is inside the islet

    for i = 1:size(Tcellsold,1)
        if(inpolygon(Tcellsold(i,1),Tcellsold(i,2),Selection.Vertices(:,1),Selection.Vertices(:,2)))
            Tcells(end+1,1) = Tcellsold(i,1);
            Tcells(end,2) = Tcellsold(i,2);
            con(end+1,1) = 0;
            if (inpolygon(Tcellsold(i,1),Tcellsold(i,2),isletselection.Vertices(:,1),isletselection.Vertices(:,2)))
                con(end,1) = 1;
            end
        end
    end

    %Retain only the alpha cells that are inside the mask

    Acells = [];

    for i = 1:size(Acellsold,1)
        if (inpolygon(Acellsold(i,1),Acellsold(i,2),Selection.Vertices(:,1),Selection.Vertices(:,2)))
            Acells(end+1,:) = Acellsold(i,:);
        end
    end

    %Separate immune cells inside and outside the islet

    TcellsOutside = [];
    TcellsInside = [];

    for i = 1:size(con)
        if con(i) == 1
            TcellsInside(end+1,:) = Tcells(i,:);
        else
            TcellsOutside(end+1,:) = Tcells(i,:);
        end
    end

    %% Construct experimental alpha-immune cell network

    %Determine number of cells
    nAcells = size(Acells,1);
    nTcells = size(Tcells,1);

    % Store list of cells with type in array
    twocells = [Acells, ones(nAcells,1); Tcells, 2*ones(nTcells,1)];

    %Define cell ranges
    Acellsrange = 1:nAcells;
    Tcellsrange = nAcells+1:nAcells+nTcells;

    ExpNetwork = zeros(size(twocells,1),size(twocells,1)); %Initialize network

    %Create network
    for i =1:size(twocells,1)
        for j = 1:size(twocells,1)
            dist = distancebetweenpoints(twocells(i,1:2),twocells(j,1:2));
            if dist<threshold && twocells(i,3) ~= twocells(j,3)
                ExpNetwork(i,j) = 1;
            end
        end
    end
    ExpNetwork = ExpNetwork - diag(diag(ExpNetwork));
    ExpNetworkgraph = graph(ExpNetwork);

    %% Compute Kavg of experimental network

    Kavgcount = []; %Define an empty array to store the number of links of each cell in an islet

    for i = Acellsrange 
        count = 0; %Initalize number of links of a cell to zero
        % Count the number of links for a cell
        for j = Tcellsrange
            if ExpNetwork(i,j) == 1
                count = count + 1;
            end
        end
        Kavgcount(end+1) = count; %Store the number of links of a cell
        Kavgcellexp(end+1,1) = count; %Store the number of links of a cell 
    end

    Kavgisletexp(end+1,1) = mean(Kavgcount); %Determine the average links of each alpha cell in an islet

    %% Create and analyze different random networks

    for k = 1:numseeds %Set seed

        %% Randomization of the positions of immune cells

        %Apply normal distribution to randomly generate a list of random distances between an immune cell and islet rim 

        Tcellsrand = []; %Initialize a variable to store randomly generated immune cells

        %For immune cells outside the islet

        if ~isempty(TcellsOutside)

            %Determine mean and standard deviation of distances of immune cells to islet

            distTOutsidetoIslet = min(pdist2(TcellsOutside,isletselection.Vertices),[],2);
            meandistOutside = mean(distTOutsidetoIslet);
            stddistOutside = std(distTOutsidetoIslet);

            %Create a normal distribution of distances

            rng(1); %Set random generator (fixed to 1)
            randdistOutside =normrnd(meandistOutside,stddistOutside,1000,1);
            randdistOutside = round(randdistOutside);

            %Randomly pick distances

            rng(k); %Set seed of random generator

            tcount = size(TcellsOutside,1);  %Determine the number of immune cells outside the islet
            indicesOutside = randi(size(randdistOutside,1),1,tcount); %Randomly pick a list of indices for distances from the normal distribution

            %Randomly assign coordinates for each random distance

            for i = 1:tcount
                dist = randdistOutside(indicesOutside(i)); %Determine the value of random distance for the selected index
                if dist > max(distCoordOutsidetoIslet)
                    dist = max(distCoordOutsidetoIslet); %Assign maximum distance if the random distance is greater than the maximum available distance 
                end
                if dist < min(distCoordOutsidetoIslet)
                    dist = min(distCoordOutsidetoIslet); %Assign minimum distance if the random distance is greater than the minimum available distance
                end
                indicesdist = find(distCoordOutsidetoIslet == dist); %Find the list of available coordinates with distance equal to the random distance

                %Randomly select a coordinate from the mask with distance equal to the random distance
                index = randi(size(indicesdist,1),1); 
                Tcellsrand(end+1,:) = coordsOutside(indicesdist(index),:); 
            end
        end

        %For immune cells inside the islet

        if ~isempty(TcellsInside)

            %Determine mean and standard deviation of distances of immune cells to islet

            distTInsidetoIslet = min(pdist2(TcellsInside,isletselection.Vertices),[],2);
            meandistInside = mean(distTInsidetoIslet);
            stddistInside = std(distTInsidetoIslet);

            %Create a normal distribution of distances

            rng(1); %Set random generator (fixed to 1)
            randdistInside =normrnd(meandistInside,stddistInside,1000,1);
            randdistInside = round(randdistInside);

            %Randomly pick distances

            rng(k); %Set seed of random generator
            tcount = size(TcellsInside,1); %Determine the number of immune cells inside the islet
            indicesInside = randi(size(randdistInside,1),1,tcount); %Randomly pick a list of indices for distances from the normal distribution

            %Randomly assign coordinates for each random distance

            for i = 1:tcount 
                dist = randdistInside(indicesInside(i)); %Determine the value of random distance for the selected index
                if dist > max(distCoordInsidetoIslet)
                    dist = max(distCoordInsidetoIslet); %Assign maximum distance if the random distance is greater than the maximum available distance
                end
                if dist < min(distCoordInsidetoIslet)
                    dist = min(distCoordInsidetoIslet); %Assign minimum distance if the random distance is greater than the minimum available distance 
                end
                indicesdist = find(distCoordInsidetoIslet == dist); %Find the list of available coordinates with distance equal to the random distance

                %Randomly select a coordinate from the mask with distance equal to the random distance
                index = randi(size(indicesdist,1),1); 
                Tcellsrand(end+1,:) = coordsInside(indicesdist(index),:); 
            end
        end

        %% Construct random alpha-immune cell network

        %Determine number of cells
        nAcells = size(Acells,1);
        nTcells = size(Tcellsrand,1);

        % Store list of cells with type in array
        twocells = [Acells, ones(nAcells,1); Tcellsrand, 2*ones(nTcells,1)];

        %Define cell ranges
        Acellsrange = 1:nAcells;
        Tcellsrange = nAcells+1:nAcells+nTcells;

        RandNetwork = zeros(size(twocells,1),size(twocells,1)); %Initialize network

        %Create network
        for i =1:size(twocells,1)
            for j = 1:size(twocells,1)
                dist = distancebetweenpoints(twocells(i,1:2),twocells(j,1:2));
                if dist<threshold && twocells(i,3) ~= twocells(j,3)
                    RandNetwork(i,j) = 1;
                end
            end
        end
        RandNetwork = RandNetwork - diag(diag(RandNetwork));
        RandNetworkgraph = graph(RandNetwork);

        %% Compute Kavg of random network

        Kavgcount = []; %Define an empty array to store the number of links of each cell in an islet

        for i = Acellsrange 
            count = 0; %Initalize number of links of a cell to zero
            % Count the number of links for a cell
            for j = Tcellsrange
                if RandNetwork(i,j) == 1
                    count = count + 1;
                end
            end
            Kavgcount(end+1) = count; %Store the number of links of a cell
            Kavgrandcell{k}(end+1,1) = count; %Store the number of links of a cell 
        end

        Kavgrandislet{k}(end+1,1) = mean(Kavgcount); %Determine the average links of each alpha cell in an islet

    end

end

%% Save the information of network analysis in one array

%Save results of network analysis of random islets

for i = 1:numseeds
    Kavgcellrand(:,i) = Kavgrandcell{i};
end
for i = 1:numseeds
    Kavgisletrand(:,i) = Kavgrandislet{i};
end

%Group files for easy export

Kavgisletwise = [Kavgisletexp mean(Kavgisletrand,2)];
Kavgcellwise = [Kavgcellexp mean(Kavgcellrand,2)];

%% Define function to calculate distance between points

function dist = distancebetweenpoints(A,B)
    dist = sqrt((A(1)-B(1))^2 + (A(2)-B(2))^2); 
end