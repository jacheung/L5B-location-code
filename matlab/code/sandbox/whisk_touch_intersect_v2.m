hilbertVar = 'angle';

selectedCells = find(cellfun(@(x) strcmp(x.meta.touchProperties.responseType,'excited'),U));

tStruct = object_location_quantification(U,selectedCells,hilbertVar,'off');

wStruct = whisking_location_quantification(U,1:numel(U),hilbertVar,'off');


if strcmp(hilbertVar,'pole')
    population_heatmap_builder(tStruct,wStruct,hilbertVar)
    
    %     saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig2\';
    %     fn = 'population_location.eps';
    %     export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
    %     fix_eps_fonts([saveDir, fn])
else
    disp('not building out population heatmaps. function not optimized for other variables')
end
%% scatter of tuning preference of whisk and touch

tUnits = cellfun(@(x) x.is_tuned==1,tStruct);
wUnits = cellfun(@(x) x.is_tuned==1,wStruct);

[~,touch_ix_idx] = intersect(find(tUnits),find(wUnits));
[~,whisk_ix_idx] = intersect(find(wUnits),find(tUnits));

touch_nonIX_idx = setdiff(1:sum(tUnits),touch_ix_idx);
whisk_nonIX_idx = setdiff(1:sum(wUnits),whisk_ix_idx);

touch_pw = cell2mat(cellfun(@(x) [x.calculations.tune_peak x.calculations.tune_left_width x.calculations.tune_right_width],tStruct(tUnits),'uniformoutput',0)') ;
whisking_pw = cell2mat(cellfun(@(x) [x.calculations.tune_peak x.calculations.tune_left_width x.calculations.tune_right_width],wStruct(wUnits),'uniformoutput',0)');

%% TOUCH TUNING CURVES
location_units = find(tUnits);
rc = numSubplots(numel(location_units));
figure(239);clf
for g = 1:numel(location_units)
    array = U{location_units(g)};
    response_window = array.meta.touchProperties.responseWindow;
    [tVar] = atTouch_sorter(array,response_window(1):response_window(2));
    
    if strcmp(hilbertVar,'phase')
        x_raw = tVar.allTouches.S_ctk(:,5);
    elseif strcmp(hilbertVar,'angle')
        x_raw = tVar.allTouches.S_ctk(:,1);
    end
    
    y_raw = mean(tVar.allTouches.R_ntk,2) .* 1000;
    
    x = x_raw(~any(isnan([x_raw y_raw]),2));
    y = y_raw(~any(isnan([x_raw y_raw]),2));
    
    f = fit(x,y,'gauss1');
    
    figure(239);subplot(rc(1),rc(2),g)
    scatter(x,y,'.k')
    yyaxis right
    hold on; plot(linspace(min(x),max(x),20),f(linspace(min(x),max(x),20)),'g');
    hold on; plot(tStruct{location_units(g)}.stim_response.values(:,1),tStruct{location_units(g)}.stim_response.values(:,2),'r')
end

%% WHISK TUNING CURVES
whisk_units = find(wUnits);
rc = numSubplots(numel(whisk_units));
num_bins = 50;
figure(239);clf
for g = 1:numel(whisk_units)
    array = U{whisk_units(g)};
    
    spikes = squeeze(array.R_ntk);
    timePostTouchToTrim = 30;
    touchOnIdx = [find(array.S_ctk(9,:,:)==1); find(array.S_ctk(12,:,:)==1)];
    touchOffIdx = [find(array.S_ctk(10,:,:)==1); find(array.S_ctk(13,:,:)==1)];
    touchOnIdx = touchOnIdx(touchOffIdx<numel(spikes)-timePostTouchToTrim+5);
    touchOffIdx = touchOffIdx(touchOffIdx<numel(spikes)-timePostTouchToTrim+5);
    touchEx_mask = ones(size(squeeze(array.S_ctk(1,:,:))));
    for i = 1:length(touchOnIdx)
        touchEx_mask(touchOnIdx(i):touchOffIdx(i)+timePostTouchToTrim) = NaN; %added 30 to add time from touch offset
    end
    touchEx_mask(1:100,:) = 1; %since bleedover from end of trials before, tihs ensure we keep end
    
    amplitude = squeeze(array.S_ctk(3,:,:));
    whisking = nan(size(squeeze(array.S_ctk(1,:,:))));
    whisking(amplitude>5)=1;
    
    whisking_mask = whisking .* touchEx_mask;
    
    
    if strcmp(hilbertVar,'phase')
        x_raw = squeeze(array.S_ctk(5,:,:)) .* whisking_mask;
    elseif strcmp(hilbertVar,'angle')
        x_raw = squeeze(array.S_ctk(1,:,:)) .* whisking_mask;
    end
    y_raw = spikes .* whisking_mask;
    
    x = x_raw(~any(isnan([x_raw(:) y_raw(:)]),2));
    y = y_raw(~any(isnan([x_raw(:) y_raw(:)]),2));
    
    if strcmp(hilbertVar,'phase')
        [sorted, sortedBy] = binslin(x,y,'equalE',num_bins+1,-pi,pi);
        y_bin = cellfun(@(x) nanmean(x),sorted) .* 1000;
        x_bin = linspace(-pi,pi,num_bins);
    elseif strcmp(hilbertVar,'angle')
        [sorted, sortedBy] = binslin(x,y,'equalE',num_bins+1,min(x),max(x));
        y_bin = cellfun(@(x) nanmean(x),sorted) .* 1000;
        x_bin = linspace(min(x),max(x),num_bins);
    end
    
    try
        f = fit(x_bin',y_bin,'gauss1');
        
        %     f = fit(x,y .* 1000,'gauss1');
        
        figure(239);subplot(rc(1),rc(2),g)
        scatter(x_bin,y_bin,'.k')
        yyaxis right
        hold on; plot(linspace(min(x),max(x),20),f(linspace(min(x),max(x),20)),'g');
        hold on; plot(wStruct{whisk_units(g)}.stim_response.values(:,1),wStruct{whisk_units(g)}.stim_response.values(:,2),'r')
        hold on; plot(wStruct{whisk_units(g)}.stim_response.bars_fit.x,wStruct{whisk_units(g)}.stim_response.bars_fit.mean,'b')
    end
end
