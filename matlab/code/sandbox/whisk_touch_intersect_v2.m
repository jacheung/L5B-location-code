%Load whisking and neural time series struct
clear
load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\excitatory_all.mat') %L5b excitatory cells recorded by Jon and Phil
% load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\interneurons.mat') %L5b inhibitory cells
%%
touchCells = find(cellfun(@(x) strcmp(x.meta.touchProperties.responseType,'excited'),U));
hilbertVar = 'phase';
tStruct = object_location_quantification(U,touchCells,hilbertVar,'off');

fileName = ['whisk_' hilbertVar '_tune'];
if exist(['C:\Users\jacheung\Dropbox\LocationCode\DataStructs\' fileName '.mat'],'file')
    load(['C:\Users\jacheung\Dropbox\LocationCode\DataStructs\' fileName '.mat'])
else
    wStruct = whisking_location_quantification(U,1:numel(U),hilbertVar,'off');
end
%% proportion of neurons that intersect bar
hilbertVar = {'angle','phase','amplitude','midpoint'};
for i = 1:numel(hilbertVar)
    tStruct = object_location_quantification(U,touchCells,hilbertVar{i},'off');
    
    fileName = ['whisk_' hilbertVar{i} '_tune'];
    if exist(['C:\Users\jacheung\Dropbox\LocationCode\DataStructs\' fileName '.mat'],'file')
        load(['C:\Users\jacheung\Dropbox\LocationCode\DataStructs\' fileName '.mat'])
    else
        wStruct = whisking_location_quantification(U,1:numel(U),hilbertVar{i},'off');
    end
    
    tUnits = cellfun(@(x) x.is_tuned==1,tStruct);
    wUnits = cellfun(@(x) x.is_tuned==1,wStruct);
    
    both_tuned_num = numel(intersect(find(tUnits),find(wUnits)));
    touch_tuned_num = numel(setdiff(1:sum(tUnits),touch_ix_idx));
    whisk_tuned_num = numel(setdiff(1:sum(wUnits),whisk_ix_idx));
    untuned_num = numel(U) - (both_tuned_num + touch_tuned_num + whisk_tuned_num);
    
    stacked_bar{i} = [ both_tuned_num touch_tuned_num whisk_tuned_num untuned_num];
    
    
end

figure(340);clf
bar(cell2mat(stacked_bar') ./ (ones(4,1).*numel(U)),'stacked')
set(gca,'xticklabel',hilbertVar,'ytick',0:.25:1)
legend({'both','touch','whisk','none'})

saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig6\';
fn = ['hilbert_proportions.eps'];
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])



%% scatter of tuning preference of whisk and touch

tUnits = cellfun(@(x) x.is_tuned==1,tStruct);
wUnits = cellfun(@(x) x.is_tuned==1,wStruct);

[~,touch_ix_idx] = intersect(find(tUnits),find(wUnits));
[~,whisk_ix_idx] = intersect(find(wUnits),find(tUnits));

touch_nonIX_idx = setdiff(1:sum(tUnits),touch_ix_idx);
whisk_nonIX_idx = setdiff(1:sum(wUnits),whisk_ix_idx);

touch_pw = cell2mat(cellfun(@(x) [x.calculations.tune_peak x.calculations.tune_left_width x.calculations.tune_right_width],tStruct(tUnits),'uniformoutput',0)') ;
whisking_pw = cell2mat(cellfun(@(x) [x.calculations.tune_peak x.calculations.tune_left_width x.calculations.tune_right_width],wStruct(wUnits),'uniformoutput',0)');

%scatter of whisking (Y) vs touch (X)
figure(3850);clf
if strcmp(hilbertVar,'pole')
    hold on; errorbar(whisking_pw(whisk_nonIX_idx,1),ones(1,length(whisk_nonIX_idx))*2.5,whisking_pw(whisk_nonIX_idx,2),whisking_pw(whisk_nonIX_idx,3),'co','horizontal')%plot only whisk tuned units
    hold on; errorbar(ones(1,length(touch_nonIX_idx))*2.5,touch_pw(touch_nonIX_idx,1),touch_pw(touch_nonIX_idx,2),touch_pw(touch_nonIX_idx,3),'bo','vertical')%plot only touch tuned units
    hold on; errorbar(whisking_pw(whisk_ix_idx,1),touch_pw(touch_ix_idx,1),touch_pw(touch_ix_idx,2),touch_pw(touch_ix_idx,3),whisking_pw(whisk_ix_idx,2),whisking_pw(whisk_ix_idx,3),'ko')

    set(gca,'xlim',[-3 3],'ylim',[-3 3],'xdir','reverse','ydir','reverse',...
        'xtick',-1:1:1,'ytick',-1:1:1)
    hold on; plot([-1 1],[-1 1],'--k')

    figure(3851);clf
    subplot(2,1,1)
    histogram(touch_pw(:,1),-3:.20:3,'facecolor','b','facealpha',1)
    set(gca,'xdir','reverse','xlim',[-3 3])
    
    subplot(2,1,2);
    histogram(whisking_pw(:,1),-3:.20:3,'facecolor','c','facealpha',1)
    set(gca,'xdir','reverse','xlim',[-3 3],'ytick',0:2:6,'ylim',[0 6])
elseif strcmp(hilbertVar,'phase')
    subplot(4,1,[1:2])
    hold on; errorbar(whisking_pw(whisk_nonIX_idx,1),ones(1,length(whisk_nonIX_idx))*-4,whisking_pw(whisk_nonIX_idx,2),whisking_pw(whisk_nonIX_idx,3),'co','horizontal','capsize',0)%plot only whisk tuned units
    hold on; errorbar(ones(1,length(touch_nonIX_idx))*-4,touch_pw(touch_nonIX_idx,1),touch_pw(touch_nonIX_idx,2),touch_pw(touch_nonIX_idx,3),'bo','vertical','capsize',0)%plot only touch tuned units
    hold on; errorbar(whisking_pw(whisk_ix_idx,1),touch_pw(touch_ix_idx,1),touch_pw(touch_ix_idx,2),touch_pw(touch_ix_idx,3),whisking_pw(whisk_ix_idx,2),whisking_pw(whisk_ix_idx,3),'ko','capsize',0)

    set(gca,'xlim',[-4 4],'ylim',[-4 4],...
        'xtick',-pi:pi:pi,'ytick',-pi:pi:pi,'xticklabel',{'-\pi',0,'\pi'},'yticklabel',{'-\pi',0,'\pi'})
    hold on; plot([-4 4],[-4 4],'--k')
    
    subplot(4,1,3)
    histogram(touch_pw(:,1),-pi:pi/4:pi,'facecolor','b','facealpha',1)
    set(gca,'xlim',[-pi pi],'ylim',[0 20])
    subplot(4,1,4);
    histogram(whisking_pw(:,1),-pi:pi/4:pi,'facecolor','c','facealpha',1)
    set(gca,'xlim',[-pi pi],'ylim',[0 20])
elseif strcmp(hilbertVar,'angle')
    subplot(4,1,[1:2])
    hold on; errorbar(whisking_pw(whisk_nonIX_idx,1),ones(1,length(whisk_nonIX_idx))*-30,whisking_pw(whisk_nonIX_idx,2),whisking_pw(whisk_nonIX_idx,3),'co','horizontal','capsize',0)%plot only whisk tuned units
    hold on; errorbar(ones(1,length(touch_nonIX_idx))*-30,touch_pw(touch_nonIX_idx,1),touch_pw(touch_nonIX_idx,2),touch_pw(touch_nonIX_idx,3),'bo','vertical','capsize',0)%plot only touch tuned units
    hold on; errorbar(whisking_pw(whisk_ix_idx,1),touch_pw(touch_ix_idx,1),touch_pw(touch_ix_idx,2),touch_pw(touch_ix_idx,3),whisking_pw(whisk_ix_idx,2),whisking_pw(whisk_ix_idx,3),'ko','capsize',0)

    set(gca,'xlim',[-40 80],'ylim',[-40 80],...
        'xtick',-40:20:80,'ytick',-40:20:80)
    hold on; plot([-40 80],[-40 80],'--k')
    
    subplot(4,1,3)
    histogram(touch_pw(:,1),-40:10:80,'facecolor','b','facealpha',1)
    set(gca,'xlim',[-40 80],'ylim',[0 20])
    subplot(4,1,4);
    histogram(whisking_pw(:,1),-40:10:80,'facecolor','c','facealpha',1)
    set(gca,'xlim',[-40 80],'ylim',[0 20])
end

figure(3850);subplot(4,1,[1 2])
lm = fitlm(touch_pw(touch_ix_idx,1),whisking_pw(whisk_ix_idx,1));
predicts = lm.predict;
[s_vals,sort_idx] = sort(touch_pw(touch_ix_idx,1));
hold on; plot(s_vals,predicts(sort_idx))

legend('whisk tuned only','touch tuned only','both tuned')
axis square
xlabel('whisk tune peak');ylabel('touch tune peak')
title(['whisk=' num2str(numel(whisk_nonIX_idx)) ', touch=' num2str(numel(touch_nonIX_idx)) ', both=' num2str(numel(touch_ix_idx))])

figure(3850);
saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig6\';
fn = [hilbertVar '_intersect.eps'];
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])


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
