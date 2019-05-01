%Load whisking and neural time series struct 
clear
load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\excitatory.mat') %L5b excitatory cells
% load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\interneurons.mat') %L5b inhibitory cells
%% Top level parameters and definitions 
touchWindow = [-25:50]; %window for analyses around touch

touchCells = touchCell(U,'off');
selectedCells = find(touchCells==1);

% Structure for quantifying tuning and evaluating decoding 
popV = touchFeatureBinned(U,touchWindow);

% Defining touch response
U = defTouchResponse(U,.99,'off');

%% Quantifying object location tuning
whichTouches = fields(popV{1});
fieldsList = fields(popV{1}.allTouches);
tuneStruct = tuningQuantification(U,popV,selectedCells,fieldsList(1),whichTouches,touchWindow,'off');

%% Quantifying whisking tuning
whisking = whisking_general(U,'off');

%% comparison 
naiveVSexpert = cellfun(@(x) strcmp(x.meta.layer,'BVL5b'),U);

wEXCOL = numel(intersect(find(tuneStruct.matrix(1,:) == 1),find(whisking.matrix(1,:) == 1))) ./ numel(find(tuneStruct.matrix(1,:) == 1)) ;
wINHOL = numel(intersect(find(tuneStruct.matrix(1,:) == 1),find(whisking.matrix(1,:) == -1))) ./ numel(find(tuneStruct.matrix(1,:) == 1));
nsOL = numel(intersect(find(tuneStruct.matrix(1,:) == 1),find(whisking.matrix(1,:) == 0))) ./ numel(find(tuneStruct.matrix(1,:) == 1));



wEXCexpert = numel(intersect(find(naiveVSexpert==1),find(whisking.matrix(1,:) == 1))) ./ sum(naiveVSexpert==1);
wEXCnaive = numel(intersect(find(naiveVSexpert==0),find(whisking.matrix(1,:) == 1))) ./ sum(naiveVSexpert==0);

wINHexpert = numel(intersect(find(naiveVSexpert==1),find(whisking.matrix(1,:) == -1))) ./ sum(naiveVSexpert==1);
wINHnaive = numel(intersect(find(naiveVSexpert==0),find(whisking.matrix(1,:) == -1))) ./ sum(naiveVSexpert==0);

nsexpert = numel(intersect(find(naiveVSexpert==1),find(whisking.matrix(1,:) == 0))) ./ sum(naiveVSexpert==1);
nsnaive = numel(intersect(find(naiveVSexpert==0),find(whisking.matrix(1,:) == 0))) ./ sum(naiveVSexpert==0);