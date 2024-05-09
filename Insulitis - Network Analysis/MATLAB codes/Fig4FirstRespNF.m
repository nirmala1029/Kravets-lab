%% This code identifies alpha-linked and non-alpha-linked beta cells, creates proximity-based networks of these cells with immune cells and computes average links of these networks
%% Credit: Nirmala Balasenthilkumaran, Kravets lab

clc
clear all
close all

%% Initialize parameters

filepath = 'D:\Nirmala\Insulitis project\All islets\'; %Update filepath 
threshold = 25; %Distance threshold for network creation
scaling = 0.4964674; %Relationship between 1 pixel and 1 um

sno = 0;

for isletno = 1:134 %List the islets to be analyzed

    sno = sno+1;

    %Load islet mask and create a shape for selection

    maskislet = logical(imread([filepath,['islet ', num2str(isletno), ' mask 0.tif']]));
    isletSelection = bwboundaries(maskislet);
    yMask = isletSelection{1}(:,1)*scaling;
    xMask = isletSelection{1}(:,2)*scaling;
    isletSelection = polyshape(xMask,yMask);

    %Load the inner boundary of islet and create a shape for selection

    maskinner = logical(imread([filepath,['islet ', num2str(isletno), ' mask -20.tif']]));
    innerSelection = bwboundaries(maskinner);
    yMask = innerSelection{1}(:,1)*scaling;
    xMask = innerSelection{1}(:,2)*scaling;
    innerSelection = polyshape(xMask,yMask);

    %Load the islet periphery mask and create a shape for selection

    maskouter = logical(imread([filepath,['islet ', num2str(isletno), ' mask 20.tif']]));
    outerSelection = bwboundaries(maskouter);
    yMask = outerSelection{1}(:,1)*scaling;
    xMask = outerSelection{1}(:,2)*scaling;
    outerSelection = polyshape(xMask,yMask);

    %Load cell coordinates of alpha, beta and immune cells

    Acellsold = table2array(readtable([filepath,['islet ', num2str(isletno), ' glucagon position.xlsx']]));
    Bcellsold = table2array(readtable([filepath,['islet ', num2str(isletno), ' insulin position.xlsx']]));
    Tcellsold = table2array(readtable([filepath,['islet ', num2str(isletno), ' f480 position.xlsx']])); %Can load other immune cells by changing the input file of this variable
    
    %Retain only those alpha and beta cells that are in the outer boundary of the islet

    Bcells = [];
    for i = 1:size(Bcellsold,1)
        if(inpolygon(Bcellsold(i,1),Bcellsold(i,2),isletSelection.Vertices(:,1),isletSelection.Vertices(:,2)))
            if(~inpolygon(Bcellsold(i,1),Bcellsold(i,2),innerSelection.Vertices(:,1),innerSelection.Vertices(:,2)))
                Bcells(end+1,1) = Bcellsold(i,1);
                Bcells(end,2) = Bcellsold(i,2);
            end
        end
    end

    Acells = [];
    for i = 1:size(Acellsold,1)
        if(inpolygon(Acellsold(i,1),Acellsold(i,2),isletSelection.Vertices(:,1),isletSelection.Vertices(:,2)))
            if(~inpolygon(Acellsold(i,1),Acellsold(i,2),innerSelection.Vertices(:,1),innerSelection.Vertices(:,2)))
                Acells(end+1,1) = Acellsold(i,1);
                Acells(end,2) = Acellsold(i,2);
            end
        end
    end

    %Retain only those immune cells that are within 20 um of islet

    Tcells = [];
    for i = 1:size(Tcellsold,1)
        if(inpolygon(Tcellsold(i,1),Tcellsold(i,2),outerSelection.Vertices(:,1),outerSelection.Vertices(:,2)))
            Tcells(end+1,1) = Tcellsold(i,1);
            Tcells(end,2) = Tcellsold(i,2);
        end
    end


    %% Create a three-cell proximity based network to identify alpha-linked beta cells

    %Determine the number of cells
    nTcells = size(Tcells,1);
    nAcells = size(Acells,1);
    nBcells = size(Bcells,1);

    % Store list of cells with type in an array
    threecells = [Tcells, ones(nTcells,1); Acells, 2*ones(nAcells,1); Bcells, 3*ones(nBcells,1)];

    % Define cell ranges
    TcellsRange = 1:nTcells;
    AcellsRange = nTcells+1:nTcells+nAcells;    
    BcellsRange = nTcells+nAcells+1:nTcells+nAcells+nBcells;

    ThreeCellsNetwork = zeros(size(threecells,1),size(threecells,1)); %Initialize network
    
    %Create network
    for i = 1:size(threecells,1)
        for j = 1:size(threecells,1)
            dist = distancebetweenpoints(threecells(i,1:2),threecells(j,1:2));
            if dist < threshold
                ThreeCellsNetwork(i,j) = 1;
            end
        end
    end
    ThreeCellsNetwork = ThreeCellsNetwork - diag(diag(ThreeCellsNetwork));
    ThreeCellsNetworkgraph = graph(ThreeCellsNetwork);

    %% Identify alpha-linked beta cells by checking if a beta cell is linked to an alpha cell

    FirstRespBeta = [];
    for i = BcellsRange
        for j = AcellsRange
            if ThreeCellsNetwork(i,j) == 1
                FirstRespBeta(end+1) = i;
                break;
            end
        end
    end

    %% Identify non-alpha linked beta cells by checking if a beta cell is linked to an alpha cell

    NonFirstRespBeta = [];

    for i = BcellsRange
        con = 0;
        for j = AcellsRange
            if ThreeCellsNetwork(i,j) == 1
                con = 1;
            end
        end
        if con == 0
            NonFirstRespBeta(end+1) = i;
        end
    end

    %% Compute Kavg of alpha-linked beta cells

    Kavgcount = []; %Define an empty array to store the number of links of each cell in an islet
    for i = FirstRespBeta
        count = 0; %Initalize number of links of a cell to zero
        for j = TcellsRange
            if ThreeCellsNetwork(i,j) == 1
                % Count the number of links for a cell
                count = count + 1; %Store the number of links of a cell
            end
        end
        Kavgcount(end+1) = count;
    end
    KavgFirstResp(sno,1) = mean(Kavgcount); %Determine the average links of each cell of a particular cell type in an islet

    %% Compute Kavg of non-alpha linked beta cells

    Kavgcount = []; %Define an empty array to store the number of links of each cell in an islet

    for i = NonFirstRespBeta
        count = 0; %Initalize number of links of a cell to zero
        for j = TcellsRange
            if ThreeCellsNetwork(i,j) == 1
                % Count the number of links for a cell
                count = count + 1; %Store the number of links of a cell
            end
        end
        Kavgcount(end+1) = count;
    end
    KavgNonFirstResp(sno,1) = mean(Kavgcount); %Determine the average links of each cell of a particular cell type in an islet

end

%% Group results for easy export

KavgBeta = [KavgFirstResp KavgNonFirstResp];

    
%% Define function to calculate distance between points

function dist = distancebetweenpoints(A,B)
    dist = sqrt((A(1)-B(1))^2 + (A(2)-B(2))^2); 
end

