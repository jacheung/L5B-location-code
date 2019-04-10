function [glmModel] = designMatrixBuilder_hilbert(selectedArray,glmnetOpt,glmModel)

basisFunction = normalize_var(normpdf(-1*glmnetOpt.bf.bfwidth:glmnetOpt.bf.bfwidth,0,glmnetOpt.bf.bfstd),0,1);

for i = 1:length(selectedArray)
    currentCell = selectedArray{i};
    
    %Defining features
    touchmat = zeros(currentCell.t,currentCell.k);
    touchIdx = [find(currentCell.S_ctk(9,:,:)==1) ;find(currentCell.S_ctk(12,:,:)==1)];
    spikes = squeeze(currentCell.R_ntk);
    phase = squeeze(currentCell.S_ctk(5,:,:));
    amplitude = squeeze(currentCell.S_ctk(3,:,:));
    midpoint = squeeze(currentCell.S_ctk(4,:,:));
    angle = squeeze(currentCell.S_ctk(1,:,:));
    curvature = squeeze(currentCell.S_ctk(17,:,:));
    touchmat(touchIdx)=1;
    
    %Get rid of touches with nan values
    nanidx = unique([find(isnan(midpoint(touchIdx))) ; find(isnan(amplitude(touchIdx))) ; find(isnan(phase(touchIdx))) ; find(isnan(curvature(touchIdx))) ] );
    touchIdx(nanidx)=[];
    
    %convoluting features with basis fucntions
    phase_conv = conv2(basisFunction,1,phase,'same');
    amplitude_conv = conv2(basisFunction,1,amplitude,'same');
    midpoint_conv = conv2(basisFunction,1,midpoint,'same');
    curvature_conv = conv2(basisFunction,1,curvature,'same');
    touchmat_conv = conv2(basisFunction,1,touchmat,'same');
    
    %INDEX BUILDER FOR LAGS
    startIdx = glmnetOpt.buildIndices'+touchIdx';
    
    %LAG AND LAGS APART
    shiftIdx = nan(length(glmnetOpt.bf.indicesToAdd ),numel(startIdx));
    touchShiftIdxRaw=nan(numel(glmnetOpt.buildIndices),length(glmnetOpt.bf.indicesToAdd ));
    for b = 1:length(glmnetOpt.bf.indicesToAdd )
        tmpIdx{b} = startIdx + glmnetOpt.bf.indicesToAdd (b);
        shiftIdx(b,:) = tmpIdx{b}(:);
        
        touchShiftIdxRaw(:,b) = (glmnetOpt.buildIndices + glmnetOpt.bf.indicesToAdd (b))==0;
        touchShiftIdxRaw(find(touchShiftIdxRaw(:,b)==1)+(-glmnetOpt.bf.bfwidth:glmnetOpt.bf.bfwidth),b) = basisFunction;
    end
    
    touchShiftIdx = touchShiftIdxRaw(:,(glmnetOpt.bf.indicesToAdd <=0));
    
    DmatX = [repmat(touchShiftIdx,length(touchIdx),1) phase_conv(shiftIdx') amplitude_conv(shiftIdx') midpoint_conv(shiftIdx') curvature_conv(shiftIdx')];
%         DmatX = [repmat(touchShiftIdx,length(touchIdx),1) phase_conv(shiftIdx') amplitude_conv(shiftIdx') midpoint_conv(shiftIdx')];
    DmatY = spikes(startIdx(:));
    
    %remove any NAN trials
    [row,~] = find(isnan(DmatX));
    trialsToRemove = unique(ceil(row./length(glmnetOpt.buildIndices)));
    trialsToRemoveIdx = (length(glmnetOpt.buildIndices)*(trialsToRemove-1)) + (1:(length(glmnetOpt.buildIndices))) ;
    DmatX(trialsToRemoveIdx,:) = [] ;
    DmatY(trialsToRemoveIdx,:) = [] ;
    
    
    rawAngle = angle(touchIdx);
    
    glmModel{i}.io.DmatX = DmatX;
    glmModel{i}.io.DmatY = DmatY; 
    glmModel{i}.io.DmatXNormalized = (DmatX - mean(DmatX)) ./ std(DmatX);
    

    glmModel{i}.basisFunctions.touch = touchShiftIdx;
    glmModel{i}.basisFunctions.features = touchShiftIdxRaw;
    
    glmModel{i}.raw.angle = rawAngle;
    glmModel{i}.raw.trimmedAngle = rawAngle;
    glmModel{i}.raw.trimmedAngle(trialsToRemove) = [];

end