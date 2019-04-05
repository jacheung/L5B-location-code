load('C:\Users\jacheung\Documents\GitHub\L5bLocationCode\matlab\datastructs\U.mat')

%% Structure for quantifying tuning 
touchWindow = [-25:50];

popV = touchFeatureBinned(U,touchWindow);
touchCells = touchCell(U);
selectedCells = find(touchCells==1);

%% Plotter for feature tuning around touch window
gaussFilt = 1;
fieldsList = fields(popV{1});
touchFeatureBinned_plotter(U,popV,selectedCells,fieldsList(1),touchWindow,gaussFilt)

%% Plotter for object location tuning
fieldsList = fields(popV{1});
tunedCellsIdx = tuningQuantification(U,popV,selectedCells,fieldsList(1),touchWindow);