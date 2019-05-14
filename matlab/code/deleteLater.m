whichTouches = fields(popV{1});
fieldsList = fields(popV{1}.allTouches);
tuneStruct = tuningQuantification(U,popV,selectedCells,fieldsList([1 3 4 5]),whichTouches,touchWindow,'off');
