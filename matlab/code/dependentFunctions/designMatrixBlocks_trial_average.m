function [glmModel] = designMatrixBlocks_trial_average(selectedArray,glmnetOpt,glmModel)

[preDecisionTouches, pd_mask] = preDecisionTouchMat(selectedArray);
varNames = {'theta','velocity','amplitude','midpoint','phase','kappa'};

for g = 1:length(selectedArray)
    array = selectedArray{g};
    varx=[1:5 17];
    
    %Input variables 
    for variableNumber=1:length(varx)
        
        [numTouches, meanFeature] = preDecisionTouchFeatures(array,preDecisionTouches{g},varx(variableNumber),glmnetOpt.touchDirection,glmnetOpt.touchOrder);
        
        glmModel{g}.io.components.(varNames{variableNumber}) = nanmean(meanFeature)';
    end
    glmModel{g}.io.components.touchNumber = numTouches';
    
    %response variables
    glmModel{g}.io.DmatY = nansum(pd_mask{g} .* squeeze(array.R_ntk))';
    
    %raw variables
    glmModel{g}.raw.pole = array.meta.motorPosition;
    
end

