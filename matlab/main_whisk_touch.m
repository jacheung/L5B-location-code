%Load whisking and neural time series struct
clear
load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\excitatory_all.mat') %L5b excitatory cells recorded by Jon and Phil
% load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\interneurons.mat') %L5b inhibitory cells
%%
touchCells = find(cellfun(@(x) strcmp(x.meta.touchProperties.responseType,'excited'),U));
hilbertVar = 'phase';
tStruct = object_location_quantification(U,touchCells,hilbertVar,'off');

fileName = ['whisk_' hilbertVar '_window'];
if exist(['C:\Users\jacheung\Dropbox\LocationCode\DataStructs\Whisking_redo\' fileName '.mat'],'file')
    load(['C:\Users\jacheung\Dropbox\LocationCode\DataStructs\Whisking_redo\' fileName '.mat'])
else
    disp('no wStruct loaded, building from scratch')
    wStruct = whisking_location_quantification(U,1:numel(U),hilbertVar,'off');
end

tUnits = cellfun(@(x) x.is_tuned==1,tStruct);
wUnits = cellfun(@(x) x.is_tuned==1,wStruct);

[~,touch_ix_idx] = intersect(find(tUnits),find(wUnits));
[~,whisk_ix_idx] = intersect(find(wUnits),find(tUnits));

touch_nonIX_idx = setdiff(1:sum(tUnits),touch_ix_idx);
whisk_nonIX_idx = setdiff(1:sum(wUnits),whisk_ix_idx);
%% proportion of neurons that intersect pie (B)
hilbertVar_list = {'angle','phase'};
stacked_bar = cell(1,numel(hilbertVar_list));
for i = 1:numel(hilbertVar_list)
    tStruct_tmp = object_location_quantification(U,touchCells,hilbertVar_list{i},'off');
    
    fileName = ['whisk_' hilbertVar_list{i} '_window'];
    if exist(['C:\Users\jacheung\Dropbox\LocationCode\DataStructs\Whisking_redo\' fileName '.mat'],'file')
        load(['C:\Users\jacheung\Dropbox\LocationCode\DataStructs\Whisking_redo\' fileName '.mat'])
    else
        wStruct = whisking_location_quantification(U,1:numel(U),hilbertVar_list{i},'off');
    end
    
    tUnits = cellfun(@(x) x.is_tuned==1,tStruct_tmp);
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
fn = ['B_intersect_proportions.eps'];
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])

%% whisk and touch tuning curves (C)
co_tuned = find(tUnits.*wUnits);
touch_only = setdiff(find(tUnits),co_tuned);
whisk_only = setdiff(find(wUnits),co_tuned);
tuned_mat = {co_tuned,touch_only,whisk_only};

num_cells_max = max(cellfun(@numel, tuned_mat));
t_tick = {'','','--'}; %tick labels for touch cells co-tune, touch only, and whisk only
w_tick = {'','--',''};
save_labels = {'co_tune','touch_only','whisk_only'};

for g = 1:numel(tuned_mat)
    %find units that are built
    built_tx = find(cellfun(@(x) isfield(x,'stim_response'),tStruct));
    built_wx = find(cellfun(@(x) isfield(x,'stim_response'),wStruct));
    tx_filt2 = find(cellfun(@(x) isfield(x.stim_response,'bars_fit'),tStruct(built_tx)));
    wx_filt2 = find(cellfun(@(x) isfield(x.stim_response,'bars_fit'),wStruct(built_wx)));

    tx_ix = intersect(built_tx(tx_filt2),tuned_mat{g});
    wx_ix = intersect(built_wx(wx_filt2),tuned_mat{g});
    
    %build xyerr values based on bars fit for touch
    tx = cellfun(@(x) x.stim_response.bars_fit.x,tStruct(tx_ix),'uniformoutput',0);
    ty = cellfun(@(x) x.stim_response.bars_fit.mean,tStruct(tx_ix),'uniformoutput',0);
    nty = cellfun(@(x) normalize_var(x,0,1),ty,'uniformoutput',0); %normalized values to [0 1]
    terr = cellfun(@(x) x.stim_response.bars_fit.confBands,tStruct(tx_ix),'uniformoutput',0);
    ty_vals = {ty,nty};
    %build xyerr values based on bars fit for whisking
    wx = cellfun(@(x) x.stim_response.bars_fit.x,wStruct(wx_ix),'uniformoutput',0);
    wy = cellfun(@(x) x.stim_response.bars_fit.mean,wStruct(wx_ix),'uniformoutput',0);
    nwy = cellfun(@(x) normalize_var(x,0,1),wy,'uniformoutput',0);
    werr = cellfun(@(x) x.stim_response.bars_fit.confBands,wStruct(wx_ix),'uniformoutput',0);
    wy_vals = {wy,nwy};
    
    %plotting in one giant elongated window
    figure(58);clf
    for k = 1:numel(wy_vals)
        for b = 1:numel(tuned_mat{g})
            if k ==1
                subplot(2,num_cells_max,b)
                try
                    shadedErrorBar(tx{b},ty_vals{k}{b},terr{b}(:,2) - ty_vals{k}{b},['r' t_tick{g}])
                end
                try
                    hold on; shadedErrorBar(wx{b},wy_vals{k}{b},werr{b}(:,2) - wy_vals{k}{b},['c' w_tick{g}])
                end
            elseif k ==2
                subplot(2,num_cells_max,b+num_cells_max)
                try
                    plot(tx{b},ty_vals{k}{b},['r' t_tick{g}])
                end
                try
                    hold on; plot(wx{b},wy_vals{k}{b},['c' w_tick{g}])
                end
            end
            if strcmp(hilbertVar,'phase')
                set(gca,'xlim',[-pi pi],'xtick',[])
                if k == 2
                    set(gca,'xtick',-pi:pi:pi,'xticklabel',{'-\pi',0,'\pi'})
                end
            end
        end
    end
    saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig6\';
    fn = [hilbertVar save_labels{g} '_curves.eps'];
    export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
    fix_eps_fonts([saveDir, fn])
end
%% Whisk x touch modulation depth (D)

t_built = cellfun(@(x) isfield(x,'stim_response'),tStruct);
w_built = cellfun(@(x) isfield(x,'stim_response'),wStruct);

co_tuned = find(tUnits.*wUnits);
co_built = find(t_built.*w_built);
touch_only = setdiff(find(tUnits),co_tuned);
whisk_only = setdiff(find(wUnits),co_tuned);

touch_mod_abs = cellfun(@(x) x.calculations.mod_idx_abs,tStruct(co_built));
whisk_mod_abs = cellfun(@(x) x.calculations.mod_idx_abs,wStruct(co_built));

touch_mod_abs_tune = cellfun(@(x) x.calculations.mod_idx_abs,tStruct(co_tuned));
whisk_mod_abs_tune = cellfun(@(x) x.calculations.mod_idx_abs,wStruct(co_tuned));

figure(8540);clf
hold on;scatter(whisk_mod_abs,touch_mod_abs,'ko')
hold on;scatter(whisk_mod_abs_tune,touch_mod_abs_tune,'filled','ko')
hold on; errorbar(mean(whisk_mod_abs_tune),mean(touch_mod_abs_tune),...
    std(touch_mod_abs_tune)./sqrt(numel(use_units)),std(touch_mod_abs_tune)./sqrt(numel(co_tuned)),...
    std(whisk_mod_abs_tune)./sqrt(numel(use_units)),std(whisk_mod_abs_tune)./sqrt(numel(co_tuned))...
    ,'ro','capsize',0)
hold on; errorbar(mean(whisk_mod_abs),mean(touch_mod_abs),...
    std(touch_mod_abs)./sqrt(numel(use_units)),std(touch_mod_abs)./sqrt(numel(co_built)),...
    std(whisk_mod_abs)./sqrt(numel(use_units)),std(whisk_mod_abs)./sqrt(numel(co_built))...
    ,'ko','capsize',0)
hold on; plot([1 100],[1 100],'--k')
set(gca,'xlim',[1 100],'ylim',[1 100],...
'yscale','log','xscale','log')
axis square

xlabel('whisk abs mod (spks/s)')
ylabel('touch abs mod (spks/s)')
axis square
[~,p,~,stats] = ttest(whisk_mod_abs,touch_mod_abs);
text(5, 5,['all built : ' 'p=' num2str(p) ', tstat=' num2str(stats.tstat) ', df=' num2str(stats.df)])

[~,p,~,stats] = ttest(whisk_mod_abs_tune,touch_mod_abs_tune);
text(5, 4,['co tuned : ' 'p=' num2str(p) ', tstat=' num2str(stats.tstat) ', df=' num2str(stats.df)])

saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig6\';
fn = [hilbertVar '_touch_x_whisk_absmod.eps'];
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])
%% Shape correlation of tuning curves (E)
intersect_correlation(tStruct,wStruct,hilbertVar)

saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig6\';
fn = [hilbertVar '_correlations.eps'];
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])

%% scatter of tuning preference of whisk and touch (F)
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
%     hold on; errorbar(whisking_pw(whisk_nonIX_idx,1),ones(1,length(whisk_nonIX_idx))*-4,whisking_pw(whisk_nonIX_idx,2),whisking_pw(whisk_nonIX_idx,3),'co','horizontal','capsize',0)%plot only whisk tuned units
%     hold on; errorbar(ones(1,length(touch_nonIX_idx))*-4,touch_pw(touch_nonIX_idx,1),touch_pw(touch_nonIX_idx,2),touch_pw(touch_nonIX_idx,3),'bo','vertical','capsize',0)%plot only touch tuned units
    hold on; errorbar(whisking_pw(whisk_ix_idx,1),touch_pw(touch_ix_idx,1),touch_pw(touch_ix_idx,2),touch_pw(touch_ix_idx,3),whisking_pw(whisk_ix_idx,2),whisking_pw(whisk_ix_idx,3),'ko','capsize',0)
    
    set(gca,'xlim',[-4 4],'ylim',[-4 4],...
        'xtick',-pi:pi:pi,'ytick',-pi:pi:pi,'xticklabel',{'-\pi',0,'\pi'},'yticklabel',{'-\pi',0,'\pi'})
    hold on; plot([-4 4],[-4 4],'--k')
    
    subplot(4,1,3)
    histogram(touch_pw(touch_nonIX_idx,1),-pi:pi/4:pi,'facecolor','b','facealpha',1)
    set(gca,'xlim',[-pi pi],'ylim',[0 20])
    subplot(4,1,4);
    histogram(whisking_pw(whisk_nonIX_idx,1),-pi:pi/4:pi,'facecolor','c','facealpha',1)
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
    histogram(touch_pw(touch_nonIX_idx,1),-30:10:80,'facecolor','b','facealpha',1,'normalization','probability')
    set(gca,'xlim',[-30 80],'ylim',[0 .5])
    subplot(4,1,4);
    histogram(whisking_pw(whisk_nonIX_idx,1),-430:10:80,'facecolor','c','facealpha',1,'normalization','probability')
    set(gca,'xlim',[-30 80],'ylim',[0 .5])
    
elseif strcmp(hilbertVar,'amplitude')
    subplot(4,1,[1:2])
    hold on; errorbar(whisking_pw(whisk_nonIX_idx,1),ones(1,length(whisk_nonIX_idx))*0,whisking_pw(whisk_nonIX_idx,2),whisking_pw(whisk_nonIX_idx,3),'co','horizontal','capsize',0)%plot only whisk tuned units
    hold on; errorbar(ones(1,length(touch_nonIX_idx))*0,touch_pw(touch_nonIX_idx,1),touch_pw(touch_nonIX_idx,2),touch_pw(touch_nonIX_idx,3),'bo','vertical','capsize',0)%plot only touch tuned units
    hold on; errorbar(whisking_pw(whisk_ix_idx,1),touch_pw(touch_ix_idx,1),touch_pw(touch_ix_idx,2),touch_pw(touch_ix_idx,3),whisking_pw(whisk_ix_idx,2),whisking_pw(whisk_ix_idx,3),'ko','capsize',0)
    
    set(gca,'xlim',[0 40],'ylim',[0 40],...
        'xtick',0:10:40,'ytick',0:10:40)
    hold on; plot([0 40],[0 40],'--k')
    
    subplot(4,1,3)
    histogram(touch_pw(:,1),0:5:40,'facecolor','b','facealpha',1,'normalization','probability')
    set(gca,'xlim',[0 40],'ylim',[0 20])
    subplot(4,1,4);
    histogram(whisking_pw(:,1),-40:10:80,'facecolor','c','facealpha',1,'normalization','probability')
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
    histogram(touch_pw(:,1),-20:10:60,'facecolor','b','facealpha',1,'normalization','probability')
    set(gca,'xlim',[-20 60],'ylim',[0 20])
    subplot(4,1,4);
    histogram(whisking_pw(:,1),-40:10:80,'facecolor','c','facealpha',1,'normalization','probability')
    set(gca,'xlim',[-20 60],'ylim',[0 20])
end

figure(3850);subplot(4,1,[1 2])
legend('whisk tuned only','touch tuned only','both tuned')
axis square
xlabel('whisk tune peak');ylabel('touch tune peak')
title(['whisk=' num2str(numel(whisk_nonIX_idx)) ', touch=' num2str(numel(touch_nonIX_idx)) ', both=' num2str(numel(touch_ix_idx))])

figure(3850);
saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig6\';
fn = [hilbertVar '_intersect.eps'];
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])


% Distance from unity line (D stats)  
t_w_vect = [touch_pw(touch_ix_idx,1) whisking_pw(whisk_ix_idx,1)];
if strcmp(hilbertVar,'angle')
    diagonal_vect = repmat([-30:80]',1,2);
elseif strcmp(hilbertVar,'phase')
    diagonal_vect = repmat((-pi:pi/20:pi)',1,2);
end

euclid_dist_all = pdist2(t_w_vect,diagonal_vect);
min_ed = min(euclid_dist_all,[],2);
[mean(min_ed) std(min_ed)]
[~,p,~,stats] = ttest(min_ed)


%% Heatmap (G)
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

%% Distance from intersect (shape)
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


