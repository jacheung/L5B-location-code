%Load whisking and neural time series struct
load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\excitatory_all.mat') %L5b excitatory cells
% load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\interneurons.mat') %L5b inhibitory cells
U = struct_editor(U);
%%
dk_corr = zeros(1,numel(U));
dt_corr = zeros(1,numel(U));
fitCoeffs = cell(1,numel(U));
fitDevExplained = zeros(1,numel(U));

for i = 1:length(U)
    
    % raw variables
    touchOn = [find(U{i}.S_ctk(9,:,:)==1)  ;find(U{i}.S_ctk(12,:,:)==1)];
    touchOff = [find(U{i}.S_ctk(10,:,:)==1)  ;find(U{i}.S_ctk(13,:,:)==1)];
    spikes = squeeze(U{i}.R_ntk(:,:,:));
    angle = squeeze(U{i}.S_ctk(1,:,:));
    touch_dkap = squeeze(U{i}.S_ctk(19,:,:));
    touch_dtheta = squeeze(U{i}.S_ctk(18,:,:));
    touch_tnum = ceil(touchOn./U{i}.t);
    poles = normalize_var(U{i}.meta.motorPosition,10,0);
    touch_poles = poles(touch_tnum);
    
    % during/at touch features (dkappa, dtheta, theta)
    dk = zeros(length(touchOn),1);
    dt = zeros(length(touchOn),1);
    for g = 1:length(touchOn)
        dk_tps = touch_dkap(touchOn(g):touchOff(g));
        [~,dk_idx] = max(abs(dk_tps));
        dk(g) = dk_tps(dk_idx);
        
        dt_tps = touch_dtheta(touchOn(g):touchOff(g));
        [~,dt_idx] = max(abs(dt_tps));
        dt(g) = dt_tps(dt_idx);
    end
    touchTheta = angle(touchOn);
    
    % touch response
    try
        response_window = U{i}.meta.touchProperties.responseWindow(1):U{i}.meta.touchProperties.responseWindow(2);
    catch
        touch_cells = cellfun(@(x) isfield(x.meta.touchProperties,'responseWindow'),U);
        rw = mean(cell2mat(cellfun(@(x) x.meta.touchProperties.responseWindow,U(touch_cells),'uniformoutput',0)'));
        response_window = round(rw(1):rw(2));
    end
    
    tresponse = sum(spikes(touchOn + response_window),2);
    
    %     figure(230);clf
    %     subplot(1,2,1);
    %     scatter(touch_poles,dk,'k.')
    %     hold on; plot([0 10],[0 0],'--k')
    %     title('max dk during touch')
    %     xlabel('normalized pole location')
    %     ylabel('max dk')
    %
    %     subplot(1,2,2);
    %     scatter(touch_poles,dt,'k.')
    %     hold on; plot([0 10],[0 0],'--k')
    %     title('max dtheta during touch')
    %     xlabel('normalized pole location')
    %     ylabel('max dt')
    
    % protraction dk/dt correlation w/ pole location
    pt = dk<0;
    dk_corr(i) = corr(touch_poles(pt)',dk(pt));
    dt_corr(i) = corr(touch_poles(pt)',dt(pt));
    
    pt_spikes = tresponse(pt);
    pt_dk = dk(pt);
    pt_theta = touchTheta(pt);
    
    norm_pt = normalize_var(pt_spikes,1,0);
    %     figure(03);clf
    %     scatter(pt_theta,pt_dk,20,repmat(norm_pt,1,3),'filled');
    %% simple glm for touch cells only
    if isfield(U{i}.meta.touchProperties,'responseWindow')
        
        glmnetOpt = glmnetSet;
        glmnetOpt.standardize = 0;
        glmnetOpt.alpha = 0.95;
        glmnetOpt.xfoldCV = 5;
        glmnetOpt.numIterations = 25;
        
        DmatX = [pt_dk pt_theta];
        norm_DmatX = (DmatX - mean(DmatX)) ./ std(DmatX);
        DmatY = pt_spikes;
        
        % test and train split
        exampleIdx = 1:length(pt_spikes);
        shuffIdx = exampleIdx(randperm(length(exampleIdx)));
        trainIdx = shuffIdx(1:round(length(pt_spikes)*.8));
        testIdx = setdiff(shuffIdx,trainIdx);
        
        trainDmatX = norm_DmatX(trainIdx',:);
        trainDmatY = DmatY(trainIdx',:);
        testDmatX = norm_DmatX(testIdx',:);
        testDmatY = DmatY(testIdx',:);
        
        % model fitting
        try
            cv = cvglmnet(trainDmatX,trainDmatY,'poisson',glmnetOpt,[],glmnetOpt.xfoldCV);
            
            % model evaluation
            fitLambda = cv.lambda_1se;
            iLambda = find(cv.lambda == fitLambda);
            fitCoeffs{i} = [cv.glmnet_fit.a0(iLambda) ; cv.glmnet_fit.beta(:,iLambda)];
            
            model = exp([ones(length(testDmatX),1),testDmatX]*fitCoeffs{i});
            mu = mean(testDmatY); % null poisson parameter
            nullLogLikelihood = sum(log(poisspdf(testDmatY,mu)));
            saturatedLogLikelihood = sum(log(poisspdf(testDmatY,testDmatY)));
            fullLogLikelihood = sum(log(poisspdf(testDmatY,model)));
            fitDevExplained(i) = (fullLogLikelihood - nullLogLikelihood)/(saturatedLogLikelihood - nullLogLikelihood);
        catch
            fitDevExplained(i) = nan;
            fitCoeffs{i} = nan(size(DmatX,2)+1,1);
        end
    else
        fitDevExplained(i) = nan;
        fitCoeffs{i} = nan(size(DmatX,2)+1,1);
    end
end



