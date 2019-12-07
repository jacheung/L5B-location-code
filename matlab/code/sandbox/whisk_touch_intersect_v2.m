%Load whisking and neural time series struct
clear
load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\excitatory_all.mat') %L5b excitatory cells recorded by Jon and Phil
% load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\interneurons.mat') %L5b inhibitory cells
%%
touchCells = find(cellfun(@(x) strcmp(x.meta.touchProperties.responseType,'excited'),U));
hilbertVar = 'angle';
tStruct = object_location_quantification(U,touchCells,hilbertVar,'off');

fileName = ['whisk_' hilbertVar '_tune'];
if exist(['C:\Users\jacheung\Dropbox\LocationCode\DataStructs\' fileName '.mat'],'file')
    load(['C:\Users\jacheung\Dropbox\LocationCode\DataStructs\' fileName '.mat'])
else
    wStruct = whisking_location_quantification(U,1:numel(U),hilbertVar,'off');
end

tUnits = cellfun(@(x) x.is_tuned==1,tStruct);
wUnits = cellfun(@(x) x.is_tuned==1,wStruct);

[~,touch_ix_idx] = intersect(find(tUnits),find(wUnits));
[~,whisk_ix_idx] = intersect(find(wUnits),find(tUnits));

touch_nonIX_idx = setdiff(1:sum(tUnits),touch_ix_idx);
whisk_nonIX_idx = setdiff(1:sum(wUnits),whisk_ix_idx);

%% correlation 
tStruct = object_location_quantification(U,1:length(U),hilbertVar,'off');
intersect_correlation(tStruct,wStruct,hilbertVar)

figure(53)
fn = [hilbertVar '_correlations.eps'];
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])

%% heatmap
population_heatmap_builder(tStruct,wStruct,hilbertVar)

saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig6\';
figure(50)
fn = [hilbertVar '_intersect_heat_raw.eps'];
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])

figure(51)
fn = [hilbertVar '_intersect_heat_squish.eps'];
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])

figure(52)
fn = [hilbertVar '_histograms.eps'];
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])

%% scatter of abs modulation depth

tuned_touch = find(tUnits);
tuned_whisk = find(wUnits); 

use_units = tuned_touch(touch_ix_idx); %same as use_units = tuned_whisk(whisk_ix_idx);

touch_mod_abs = cellfun(@(x) x.calculations.mod_idx_abs,tStruct(use_units));
whisk_mod_abs = cellfun(@(x) x.calculations.mod_idx_abs,wStruct(use_units));
max_val = ceil(max([whisk_mod_abs touch_mod_abs])./10) .* 10;


figure(8540);clf
scatter(whisk_mod_abs,touch_mod_abs,'ko')
 hold on; errorbar(mean(whisk_mod_abs),mean(touch_mod_abs),...
     std(touch_mod_abs)./sqrt(sum(tuned_units)),std(touch_mod_abs)./sqrt(sum(tuned_units)),...
     std(whisk_mod_abs)./sqrt(sum(tuned_units)),std(whisk_mod_abs)./sqrt(sum(tuned_units))...
     ,'ro','capsize',0)
hold on; plot([0 max_val],[0 max_val],'--k')
set(gca,'xlim',[0 max_val],'ylim',[0 max_val],'xtick',0:25:100,'ytick',0:25:100)
axis square

xlabel('whisk abs mod (spks/s)')
ylabel('touch abs mod (spks/s)')
axis square
[~,p,~,stats] = ttest(whisk_mod_abs,touch_mod_abs);
title(['p=' num2str(p) ', tstat=' num2str(stats.tstat) ', df=' num2str(stats.df)])

saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig5\';
fn = 'touch_x_whisk_absmod.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])

%% proportion of neurons that intersect pie
hilbertVar_list = {'angle','phase','amplitude','midpoint'};
stacked_bar = cell(1,numel(hilbertVar_list));
for i = 1:numel(hilbertVar_list)
    tStruct = object_location_quantification(U,touchCells,hilbertVar_list{i},'off');
    
    fileName = ['whisk_' hilbertVar_list{i} '_tune'];
    if exist(['C:\Users\jacheung\Dropbox\LocationCode\DataStructs\' fileName '.mat'],'file')
        load(['C:\Users\jacheung\Dropbox\LocationCode\DataStructs\' fileName '.mat'])
    else
        wStruct = whisking_location_quantification(U,1:numel(U),hilbertVar_list{i},'off');
    end
    
    tUnits = cellfun(@(x) x.is_tuned==1,tStruct);
    wUnits = cellfun(@(x) x.is_tuned==1,wStruct);
    
    [~,touch_ix_idx] = intersect(find(tUnits),find(wUnits));
    [~,whisk_ix_idx] = intersect(find(wUnits),find(tUnits));
    
    both_tuned_num = numel(intersect(find(tUnits),find(wUnits)));
    touch_tuned_num = numel(setdiff(1:sum(tUnits),touch_ix_idx));
    whisk_tuned_num = numel(setdiff(1:sum(wUnits),whisk_ix_idx));
    untuned_num = numel(U) - (both_tuned_num + touch_tuned_num + whisk_tuned_num);
    
    stacked_bar{i} = [ both_tuned_num touch_tuned_num whisk_tuned_num untuned_num];
    
end

figure(340);clf
for b = 1:numel(stacked_bar)
    subplot(2,2,b)
    pie(stacked_bar{b})
    title(hilbertVar_list{b})
end
legend({'both','touch','whisk','none'})

saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig6\';
fn = ['hilbert_proportions.eps'];
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])


%% scatter of tuning preference of whisk and touch
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
%     hold on; errorbar(whisking_pw(whisk_nonIX_idx,1),ones(1,length(whisk_nonIX_idx))*-30,whisking_pw(whisk_nonIX_idx,2),whisking_pw(whisk_nonIX_idx,3),'co','horizontal','capsize',0)%plot only whisk tuned units
%     hold on; errorbar(ones(1,length(touch_nonIX_idx))*-30,touch_pw(touch_nonIX_idx,1),touch_pw(touch_nonIX_idx,2),touch_pw(touch_nonIX_idx,3),'bo','vertical','capsize',0)%plot only touch tuned units
    hold on; errorbar(whisking_pw(whisk_ix_idx,1),touch_pw(touch_ix_idx,1),touch_pw(touch_ix_idx,2),touch_pw(touch_ix_idx,3),whisking_pw(whisk_ix_idx,2),whisking_pw(whisk_ix_idx,3),'ko','capsize',0)
    
    set(gca,'xlim',[-30 80],'ylim',[-30 80],...
        'xtick',-40:20:80,'ytick',-40:20:80)
    hold on; plot([-40 80],[-40 80],'--k')
    
    subplot(4,1,3)
    histogram(touch_pw(touch_nonIX_idx,1),-30:10:80,'facecolor','b','facealpha',1)
    set(gca,'xlim',[-30 80],'ylim',[0 20])
    subplot(4,1,4);
    histogram(whisking_pw(whisk_nonIX_idx,1),-430:10:80,'facecolor','c','facealpha',1)
    set(gca,'xlim',[-30 80],'ylim',[0 20])
    
elseif strcmp(hilbertVar,'amplitude')
    subplot(4,1,[1:2])
    hold on; errorbar(whisking_pw(whisk_nonIX_idx,1),ones(1,length(whisk_nonIX_idx))*0,whisking_pw(whisk_nonIX_idx,2),whisking_pw(whisk_nonIX_idx,3),'co','horizontal','capsize',0)%plot only whisk tuned units
    hold on; errorbar(ones(1,length(touch_nonIX_idx))*0,touch_pw(touch_nonIX_idx,1),touch_pw(touch_nonIX_idx,2),touch_pw(touch_nonIX_idx,3),'bo','vertical','capsize',0)%plot only touch tuned units
    hold on; errorbar(whisking_pw(whisk_ix_idx,1),touch_pw(touch_ix_idx,1),touch_pw(touch_ix_idx,2),touch_pw(touch_ix_idx,3),whisking_pw(whisk_ix_idx,2),whisking_pw(whisk_ix_idx,3),'ko','capsize',0)
    
    set(gca,'xlim',[0 40],'ylim',[0 40],...
        'xtick',0:10:40,'ytick',0:10:40)
    hold on; plot([0 40],[0 40],'--k')
    
    subplot(4,1,3)
    histogram(touch_pw(:,1),0:5:40,'facecolor','b','facealpha',1)
    set(gca,'xlim',[0 40],'ylim',[0 20])
    subplot(4,1,4);
    histogram(whisking_pw(:,1),-40:10:80,'facecolor','c','facealpha',1)
    set(gca,'xlim',[0 40],'ylim',[0 20])
elseif strcmp(hilbertVar,'midpoint')
    subplot(4,1,[1:2])
    hold on; errorbar(whisking_pw(whisk_nonIX_idx,1),ones(1,length(whisk_nonIX_idx))*-20,whisking_pw(whisk_nonIX_idx,2),whisking_pw(whisk_nonIX_idx,3),'co','horizontal','capsize',0)%plot only whisk tuned units
    hold on; errorbar(ones(1,length(touch_nonIX_idx))*-20,touch_pw(touch_nonIX_idx,1),touch_pw(touch_nonIX_idx,2),touch_pw(touch_nonIX_idx,3),'bo','vertical','capsize',0)%plot only touch tuned units
    hold on; errorbar(whisking_pw(whisk_ix_idx,1),touch_pw(touch_ix_idx,1),touch_pw(touch_ix_idx,2),touch_pw(touch_ix_idx,3),whisking_pw(whisk_ix_idx,2),whisking_pw(whisk_ix_idx,3),'ko','capsize',0)
    
    set(gca,'xlim',[-20 60],'ylim',[-20 60],...
        'xtick',-20:20:60,'ytick',-20:20:60)
    hold on; plot([-20 60],[-20 60],'--k')
    
    subplot(4,1,3)
    histogram(touch_pw(:,1),-20:10:60,'facecolor','b','facealpha',1)
    set(gca,'xlim',[-20 60],'ylim',[0 20])
    subplot(4,1,4);
    histogram(whisking_pw(:,1),-40:10:80,'facecolor','c','facealpha',1)
    set(gca,'xlim',[-20 60],'ylim',[0 20])
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
saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig5\';
fn = [hilbertVar '_intersect.eps'];
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])


%% Distance from intersect
tUnits = find(cellfun(@(x) x.is_tuned==1,tStruct));
wUnits = find(cellfun(@(x) x.is_tuned==1,wStruct));

[intersected_units,touch_ix_idx] = intersect(tUnits,wUnits);
[intersected_whisk_units,whisk_ix_idx] = intersect(wUnits,tUnits);

touch_nonIX_idx = setdiff(1:sum(tUnits),touch_ix_idx);
whisk_nonIX_idx = setdiff(1:sum(wUnits),whisk_ix_idx);

peak_response = cellfun(@(x) x.calculations.tune_peak,wStruct(intersected_whisk_units));

figure(30);clf
interp_norm_y = cell(1,length(peak_response));
for g = 1:length(intersected_units)
    sr = tStruct{intersected_units(g)}.stim_response;
    centered_x = sr.values(:,1) - peak_response(g) ;
    norm_y = norm_new(sr.values(:,2));
    
    if strcmp(hilbertVar,'angle')
        interp_centered_x = -60:1:60;
    elseif strcmp(hilbertVar,'phase')
        interp_centered_x = linspace(-pi,pi,13);
    end
    raw_x = round(centered_x,2);
    [~,idx] = unique(raw_x);
    interp_norm_y{g} = interp1(raw_x(idx),norm_y(idx),interp_centered_x);
    hold on; plot(interp_centered_x,interp_norm_y{g},'color',[.8 .8 .8])
end

pop_mean = nanmean(cell2mat(interp_norm_y'));
pop_sem = nanstd(cell2mat(interp_norm_y')) ./ sqrt(sum(~isnan(cell2mat(interp_norm_y'))));
hold on; shadedErrorBar(interp_centered_x,pop_mean,pop_sem,'k')

if strcmp(hilbertVar,'angle')
    set(gca,'xlim',[-60 60],'ytick',0:.5:1,'xtick',-60:20:60)
elseif strcmp(hilbertVar,'phase')
    set(gca,'xlim',[-pi pi],'ytick',0:.5:1,'xtick',-pi:pi:pi,'xticklabel',{'-\pi',0,'\pi'})
elseif strcmp(hilbertVar,'pole')
    set(gca,'xlim',[-2 2],'xdir','reverse','ytick',0:.5:1,'xtick',-2:1:2)
end

xlabel('distance relative to whisk peak')
axis square
saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig6\';
fn = [hilbertVar '_distance_from_whisk.eps'];
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])

%% whisk and touch tuning curves
figure(8);clf
figure(9);clf
figure(10);clf

rc = numSubplots(numel(U));

for rec = 1:length(U)
    [tVar] = atTouch_sorter(U{rec},-25:50);
    touch_response = nanmean(tVar.allTouches.R_ntk) * 1000;
    
    %     figure(8)
    %     subplot(rc(1),rc(2),rec)
    %     bar(-25:50,touch_response,'k')
    %
    try
        whisk_xy = [wStruct{rec}.stim_response.bars_fit.x' wStruct{rec}.stim_response.bars_fit.mean...
            wStruct{rec}.stim_response.bars_fit.confBands - wStruct{rec}.stim_response.bars_fit.mean];
        norm_whisk =[wStruct{rec}.stim_response.bars_fit.x' normalize_var(whisk_xy(:,2),0,1)];
    end
    
    if wStruct{rec}.is_tuned==1
%         figure(9)
%         subplot(rc(1),rc(2),rec)
%         hold on; shadedErrorBar(whisk_xy(:,1),whisk_xy(:,2),whisk_xy(:,3),'c')
        figure(10);
        subplot(rc(1),rc(2),rec)
        hold on; plot(norm_whisk(:,1),norm_whisk(:,2),'c')
    else
%         figure(9);
%         subplot(rc(1),rc(2),rec)
%         hold on; shadedErrorBar(whisk_xy(:,1),whisk_xy(:,2),whisk_xy(:,3),'c--')
        figure(10);
        subplot(rc(1),rc(2),rec)
        hold on; plot(norm_whisk(:,1),norm_whisk(:,2),'c--')
    end
    
    %if touch cell...
    if strcmp(U{rec}.meta.touchProperties.responseType,'excited')
        responseWindow = U{rec}.meta.touchProperties.responseWindow;
        peak_y = ceil(max(touch_response)/10)*10;
        
        loc_xy = [tStruct{rec}.stim_response.bars_fit.x' tStruct{rec}.stim_response.bars_fit.mean...
            tStruct{rec}.stim_response.bars_fit.confBands - tStruct{rec}.stim_response.bars_fit.mean];
        norm_loc = [loc_xy(:,1) normalize_var(loc_xy(:,2),0,1)];
        
        %         figure(8)
        %         subplot(rc(1),rc(2),rec)
        %         hold on; plot(responseWindow,[peak_y peak_y],'r')
        %         set(gca,'ylim',[0 peak_y * 1.2])
        
        
        if tStruct{rec}.is_tuned==1
%             figure(9);
%             subplot(rc(1),rc(2),rec)
%             shadedErrorBar(loc_xy(:,1),loc_xy(:,2),loc_xy(:,3),'b')
            figure(10);
            subplot(rc(1),rc(2),rec)
            plot(norm_loc(:,1),norm_loc(:,2),'b')
        else
%             figure(9);
%             subplot(rc(1),rc(2),rec)
%             shadedErrorBar(loc_xy(:,1),loc_xy(:,2),loc_xy(:,3),'b--')
            figure(10);
            subplot(rc(1),rc(2),rec)
            plot(norm_loc(:,1),norm_loc(:,2),'b--')
        end
        
    end
    
end

% figure(9);
% saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig6\';
% fn = 'whisk_touch_tuning_curves.eps';
% export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
% fix_eps_fonts([saveDir, fn])

figure(10);
saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig6\';
fn = 'whisk_touch_tuning_curves_normalized.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])

