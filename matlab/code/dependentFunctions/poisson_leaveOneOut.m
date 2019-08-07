function glmModel = poisson_leaveOneOut(glmModel,glmnetOpt)

for i=1:length(glmModel)
    coefficients = glmModel{i}.coeffs.raw;
    numVariables = 1:size(coefficients,1);
    
    fitDevExplained = nan(size(coefficients,2),glmnetOpt.numIterations); 
    
    for g = 1:size(coefficients,2)
        sel_coefficients = coefficients(setdiff(numVariables,g),:);
        
        for d = 1:glmnetOpt.numIterations
            coeff_check = sel_coefficients(:,d);
            DmatX_check = mdl.predicted.inputX{d}(:,setdiff(numVariables,g));
            
            model = exp([ones(length(DmatX_check),1),DmatX_check]*coeff_check);
            mu = mean(testDmatY); % null poisson parameter
            nullLogLikelihood = sum(log(poisspdf(testDmatY,mu)));
            saturatedLogLikelihood = sum(log(poisspdf(testDmatY,testDmatY)));
            fullLogLikelihood = sum(log(poisspdf(testDmatY,model)));
            fitDevExplained(g,d) = (fullLogLikelihood - nullLogLikelihood)/(saturatedLogLikelihood - nullLogLikelihood);
            %devianceFullNull = 2*(fullLogLikelihood - nullLogLikelihood);
        end
    end
    
    glmModel{i}.coeff.leaveOneOutNames = {'row = feature','column = numIteration'};
    glmModel{i}.coeff.leaveOneOut = fitDevExplained;
    
end

