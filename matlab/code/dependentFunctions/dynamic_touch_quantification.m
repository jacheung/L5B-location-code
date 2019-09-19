function [tuneStruct] = dynamic_touch_quantification(uberarray,selectedCells,dynamic_touch_variable,displayOpt)

if (nargin < 4), displayOpt = 'on'; end
willdisplay = ~(strcmp(displayOpt,'nodisplay') | strcmp(displayOpt,'n') ...
    | strcmp(displayOpt,'off'));

if willdisplay
    rc = numSubplots(numel(selectedCells));
    figure(25);clf
end

dt_var_list = {'dkappa','dtheta'};
if ~any(contains(dt_var_list,dynamic_touch_variable))
    error('use dynamic touch variables of "dkappa" or "dtheta"');
end

numTouchesPerBin = 75; %number of touches to assign in each bin for quantification.
alpha_value = .05; %p-value threshold to determine whether a cell is OL tuned or not
smoothing_param = 5; %smoothing parameter for smooth f(x) in shadedErrorBar
min_bins = 5; %minimum number of angle bins to consider quantifying

tuneStruct = cell(1,numel(uberarray)); %dynamic touch tuning structure
dt = cell(1,numel(uberarray)); %dynamic touch variable structure
%% BUILD dynamic touch data structure

for rec = 1:length(selectedCells)
    array = uberarray{selectedCells(rec)};
    spikes = squeeze(array.R_ntk);
    touchOnIdx = [find(array.S_ctk(9,:,:)==1); find(array.S_ctk(12,:,:)==1)];
    touchOffIdx = [find(array.S_ctk(10,:,:)==1); find(array.S_ctk(13,:,:)==1)];
    touch_response_window = array.meta.touchProperties.responseWindow(1): array.meta.touchProperties.responseWindow(2);
    
    for i = 1:length(touchOnIdx)
        touch_window_idx = touchOnIdx(i):touchOffIdx(i);
        kwin=array.S_ctk(6,touch_window_idx); %get values in touch window
        [~ ,maxidx] = max(abs(kwin)); %find idx of max kappa within each touch window, neg or pos
        dt{selectedCells(rec)}.stim.dkappa(i)=kwin(maxidx); %use idx to pull max kappa
        dk_response_window = touch_window_idx(maxidx) + touch_response_window ;
        dt{selectedCells(rec)}.response.dkappa(i,:) = spikes(dk_response_window);
        
        dtheta=array.S_ctk(18,touchOnIdx(i):touchOffIdx(i)); %get values in touch window
        [~ ,maxidx] = max(abs(dtheta)); %find idx of max theta within each touch window, neg or pos
        dt{selectedCells(rec)}.stim.dtheta(i)=dtheta(maxidx); %use idx to pull max theta
        dt_response_window = touch_window_idx(maxidx) + touch_response_window ;
        dt{selectedCells(rec)}.response.dtheta(i,:) = spikes(dt_response_window);
    end
end

%% QUANTIFY dynamic touch data structure

for rec = 1:length(selectedCells)
    
    curr_dt = dt{selectedCells(rec)};
    numBins = round(size(curr_dt.stim.(dynamic_touch_variable),2) ./ numTouchesPerBin);
    
    [sorted, sortedBy] = binslin(curr_dt.stim.(dynamic_touch_variable),mean(curr_dt.response.(dynamic_touch_variable),2) * 1000,'equalN',numBins);
    
    if numBins >= 2
        x = cellfun(@nanmedian,sortedBy);
        y = smooth(cellfun(@nanmean,sorted),smoothing_param);
        err = smooth(cellfun(@(x)nanstd(x)./sqrt(numel(x)),sorted),smoothing_param);
        if willdisplay
            subplot(rc(1),rc(2),rec)
            shadedErrorBar(x,y,err)
            title(num2str((max(y) - min(y)) ./ (max(y) + min(y))))
        end
        
        tuneStruct{selectedCells(rec)}.calculations.mod_idx_relative = (max(y) - min(y)) ./ (max(y) + min(y));
    else
        tuneStruct{selectedCells(rec)}.calculations.mod_idx_relative = 0;
    end
    
    tuneStruct{selectedCells(rec)}.dynamic_touch_raw = curr_dt;
end


