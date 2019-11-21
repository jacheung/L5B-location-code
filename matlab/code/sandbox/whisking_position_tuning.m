saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig4\';
%% example whisker trace

for k = 83
    currentCell = U{k};
    
%      trialSample = [datasample(1:currentCell.k,5)];
     trialSample = [134 136 datasample(1:currentCell.k,1)]; %CELL 83
%         trialSample = [99 136 datasample(1:currentCell.k,3)]; %% CELL 22
%         trialSample = [119 98 43 63 datasample(1:currentCell.k,1)]; %% CELL 33 EXAMPLE 119 SO GOOD!CELL TUNED TO ANGLE 10 . Dam, it's also phase tuned. LULz
%         trialSample = [60 1 datasample(1:currentCell.k,1)]; %% CELL 74 maybe 96 too
    figure(2310);clf
    for g = 1:length(trialSample)
        
        whisk = currentCell.S_ctk(1,:,trialSample(g));
        midpoint = currentCell.S_ctk(4,:,trialSample(g)); 
        amp = currentCell.S_ctk(3,:,trialSample(g));
        amp_whisk = whisk;
        amp_whisk(amp<5) = nan;
        touchesOn = [find(currentCell.S_ctk(9,:,trialSample(g))==1) find(currentCell.S_ctk(12,:,trialSample(g))==1)];
        touchesOff = [find(currentCell.S_ctk(10,:,trialSample(g))==1) find(currentCell.S_ctk(13,:,trialSample(g))==1)];

        if ~isempty(touchesOn)
            for p = 1:numel(touchesOn)
                amp_whisk(touchesOn(p):touchesOff(p)+30) = nan;
            end
        end
        
        spikes = find(currentCell.R_ntk(1,:,trialSample(g))==1);
        subplot(length(trialSample),1,g)
        plot(1:currentCell.t,whisk,'k')
        hold on;plot(1:currentCell.t,amp_whisk,'c')
        
        hold on; scatter(spikes,whisk(spikes),'ko','filled','markerfacecolor',[.8 .8 .8])
        hold on; scatter(spikes,amp_whisk(spikes),'ko','filled')

        if ~isempty(touchesOn) && ~isempty(spikes) 
            for f = 1:length(touchesOn)
                tp = touchesOn(f):touchesOff(f);
                hold on; plot(tp,whisk(tp),'r');
            end
        end
        
        title(['cell num ' num2str(k) ' at trial num ' num2str(trialSample(g))])
        set(gca,'ylim',[-10 50])
    end
    
end

        fn = 'example_traces .eps';
        export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
        fix_eps_fonts([saveDir, fn])
%% whisk x quiet
masks = cellfun(@(x) maskBuilder(x),U,'uniformoutput',0);
whisking_spks_mat = cellfun(@(x,y) squeeze(x.R_ntk).*y.whisking .*y.touch,U,masks,'uniformoutput',0);
quiet_spks_mat = cellfun(@(x,y) squeeze(x.R_ntk).*y.quiet .*y.touch,U,masks,'uniformoutput',0);
whisking_tp = cellfun(@(x,y) nansum(nansum(y.whisking .*y.touch)),U,masks);
quiet_tp = cellfun(@(x,y) nansum(nansum(y.quiet .*y.touch)),U,masks);

[~,p] = cellfun(@(x,y) ttest2(x(:),y(:)),whisking_spks_mat,quiet_spks_mat);
fr_whisk = (cellfun(@(x) nansum(x(:)),whisking_spks_mat)./whisking_tp)*1000;
fr_quiet = (cellfun(@(x) nansum(x(:)),quiet_spks_mat)./quiet_tp)*1000;

red_dots = intersect(find(p<.05),find(fr_whisk>fr_quiet));
blue_dots = intersect(find(p<.05),find(fr_whisk<fr_quiet));
gray_dots = setdiff(1:numel(masks),[red_dots blue_dots]);

figure(480);clf
% scatter(fr_quiet(red_dots),fr_whisk(red_dots),'filled','r')
% hold on; scatter(fr_quiet(blue_dots),fr_whisk(blue_dots),'filled','b')
% hold on; scatter(fr_quiet(gray_dots),fr_whisk(gray_dots),'filled','markerfacecolor',[.8 .8 .8])
scatter(fr_quiet(p<.05),fr_whisk(p<.05),'filled','k')
hold on; scatter(fr_quiet(gray_dots),fr_whisk(gray_dots),'ko')
axis square
hold on; plot([0.01 100],[0.01 100],'--k')
set(gca,'xlim',[0 100],'ylim',[0 100],'xscale','log','yscale','log')
xlabel('quiet FR');ylabel('whisking FR')
title(['red=' num2str(numel(red_dots)) ' blue=' num2str(numel(blue_dots)) ' n.s.=' num2str(numel(gray_dots))])

fn = 'whisk_quiet_unity.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])


%% tuned units
hilbert_feature = 'angle';
fileName = ['whisk_' hilbert_feature '_tune'];
if exist(['C:\Users\jacheung\Dropbox\LocationCode\DataStructs\' fileName '.mat'],'file')
    load(['C:\Users\jacheung\Dropbox\LocationCode\DataStructs\' fileName '.mat'])
    angle_whisk = wStruct;
else
    angle_whisk = whisking_location_quantification(U,1:numel(U),hilbertVar,'off');
end
tuned_units = cellfun(@(x) x.is_tuned,angle_whisk)==1;
% 
% fn = 'whisker_position_indiv.eps';
% export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
% fix_eps_fonts([saveDir, fn])

%% scatter of absolute modulation index x FR
built = find(cellfun(@(x) isfield(x,'stim_response'),angle_whisk)); 

masks = cellfun(@(x) maskBuilder(x),U(built),'uniformoutput',0);

whisking_spks_mat = cellfun(@(x,y) squeeze(x.R_ntk).*y.whisking .*y.touch,U(built),masks,'uniformoutput',0);
whisking_tp = cellfun(@(x,y) nansum(nansum(y.whisking .*y.touch)),U(built),masks);
whisk_fr = (cellfun(@(x) nansum(x(:)),whisking_spks_mat)./whisking_tp)*1000;
mod_idx_all_units = cellfun(@(x) x.calculations.mod_idx_abs,angle_whisk(built));
mod_idx_relative_all_units = cellfun(@(x) x.calculations.mod_idx_relative,angle_whisk(built));
% tuned_units = cellfun(@(x) x.is_tuned,angle_whisk(built))==1;
tuned_units = cellfun(@(x) x.is_tuned==1,angle_whisk);
non_tuned = ~tuned_units;


figure(98);clf
subplot(4,2,[1 3])
scatter(whisk_fr(non_tuned),mod_idx_all_units(non_tuned),'ko')
hold on; scatter(whisk_fr(tuned_units),mod_idx_all_units(tuned_units),'filled','ko')
hold on; plot([0.1 100],[0.1 100],'--k')
set(gca,'xscale','log','yscale','log','ylim',[0.1 100],'xlim',[0.1 100])
xlabel('whisking firing rate');
ylabel('absolute modulation index')
axis square

subplot(4,2,[2 4]);
prop_tuned = sum(tuned_units)./ numel(tuned_units); 
prop_untuned = sum(non_tuned)./numel(tuned_units); 
pie([prop_tuned prop_untuned])

subplot(4,2,[5 6]) 
histogram(mod_idx_all_units,[0.1 .25 .5 1 2.5 5 10 25 50 100],'FaceColor','k','FaceAlpha',1)
title('mod idx')
set(gca,'xscale','log','xlim',[.1 100])

subplot(4,2,[7 8])
hold on; histogram(whisk_fr,[0.1 .25 .5 1 2.5 5 10 25 50 100],'FaceColor','k','FaceAlpha',1)
title('whisk fr')
set(gca,'xscale','log','xlim',[.1 100])

fn = 'fr_x_modulation_index.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])

%% CDF of sampling
built = find(cellfun(@(x) isfield(x,'stim_response'),angle_whisk)); 
bars_built = cellfun(@(x) isfield(x.stim_response,'bars_fit'),angle_whisk(built));
full_build = built(bars_built);

masks = cellfun(@(x) maskBuilder(x),U(full_build),'uniformoutput',0);

thetas = cellfun(@(x,y) squeeze(x.S_ctk(1,:,:)) .* y.whisking,U(full_build),masks,'uniformoutput',0);

hist_thetas = cellfun(@(x) histcounts(x(:),-25.5:80.5),thetas,'uniformoutput',0);

cdf_sums = cell2mat(cellfun(@(x) cumsum(x),hist_thetas','uniformoutput',0))';
cdf_mat = cdf_sums./(cdf_sums(end,:));

mean_time_per_trial_whisking = cellfun(@(x) nansum(x.whisking(:))./size(x.whisking,2)./1000,masks); % in seconds
num_bins_per_cell = cellfun(@(x) numel(x.stim_response.raw_stim),angle_whisk(full_build)); 

figure(23);clf
subplot(2,3,[1 2 4 5])
plot(repmat(-25:80,size(cdf_mat,2),1)',cdf_mat,'color',[.9 .9 .9])
hold on; shadedErrorBar(-25:80,nanmean(cdf_mat,2),nanstd(cdf_mat,[],2)./sqrt(numel(full_build)))
axis square
set(gca,'ytick',0:.2:1,'xtick',-20:20:80)
xlabel('whisker angle during whisking')
ylabel('cumulative proportion of whisks')

subplot(2,3,3);
scatter(ones(1,numel(full_build)),mean_time_per_trial_whisking,'k','markeredgecolor',[.9 .9 .9])
hold on; errorbar(1,mean(mean_time_per_trial_whisking),std(mean_time_per_trial_whisking),'ko')
set(gca,'ylim',[0 3],'ytick',0:1:3,'xtick',[])
ylabel('time spent whisking per trial (s)');
axis square
title(['mean=' num2str(mean(mean_time_per_trial_whisking)) ' , SD=' num2str(std(mean_time_per_trial_whisking))])

subplot(2,3,6);
scatter(ones(1,numel(full_build)),num_bins_per_cell,'k','markeredgecolor',[.9 .9 .9])
hold on; errorbar(1,mean(num_bins_per_cell),std(num_bins_per_cell),'ko')
set(gca,'ylim',[0 80],'ytick',0:25:100,'xtick',[])
ylabel('number of whisking bins');
axis square
title(['mean=' num2str(mean(num_bins_per_cell)) ' , SD=' num2str(std(num_bins_per_cell))])

fn = 'whisk_sampling.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])

%% phase map
hilbert_feature = 'phase'
fileName = ['whisk_' hilbert_feature '_tune'];
if exist(['C:\Users\jacheung\Dropbox\LocationCode\DataStructs\' fileName '.mat'],'file')
    load(['C:\Users\jacheung\Dropbox\LocationCode\DataStructs\' fileName '.mat'])
    phase_whisk = wStruct;
else
    phase_whisk = whisking_location_quantification(U,1:numel(U),hilbertVar,'off');
end

tuned_units = cellfun(@(x) x.is_tuned,angle_whisk)==1;

preferred_phase = cellfun(@(x) x.calculations.tune_peak,phase_whisk(tuned_units));
phase_mod = cellfun(@(x) x.calculations.mod_idx_relative,phase_whisk(tuned_units)); 
retractions = sign(preferred_phase)<0;
conversions = zeros(1,numel(preferred_phase));
conversions(retractions)=180;
deg_per_pie_unit =  180/pi; 

preferred_degs = (abs(preferred_phase) .* deg_per_pie_unit) + conversions;

figure(380);clf
polarplot(preferred_phase,phase_mod,'o');

saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig4\';
fn = 'phase_map.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])
%% pop map 
whisk_angle_tuned= cellfun(@(x) x.is_tuned==1,angle_whisk);

bars_io = cellfun(@(x) [x.stim_response.bars_stim' x.stim_response.bars_fit.mean ],angle_whisk(whisk_angle_tuned),'uniformoutput',0);

whisk_sampled = round([min(cellfun(@(x) min(x(:,1)),bars_io)):max(cellfun(@(x) max(x(:,1)),bars_io))]);

interped_resp = cellfun(@(x) interp1(x(:,1),x(:,2),whisk_sampled),bars_io,'uniformoutput',0);
% norm_resp = cellfun(@(x) normalize_var(x,0,1),interped_resp,'uniformoutput',0);
norm_resp = cellfun(@(x) norm_new(x),interped_resp,'uniformoutput',0);

[~,maxidx] = cellfun(@(x) max(x),norm_resp);
[~,sorted_idx] = sort(maxidx); 

sort_whisk_tuning = cell2mat(norm_resp(sorted_idx)'); 

figure(45);clf
pcolor(flipud(cell2mat(norm_resp(sorted_idx)'))); %flipud because pcolor flips compared to imagesc
% set(gca,'xtick',1:10:numel(whisk_sampled),'xticklabel',whisk_sampled(1:10:end))
set(gca,'xtick',4:10:numel(whisk_sampled),'xticklabel',whisk_sampled(4:10:end))
colormap turbo
colorbar
caxis([0 1])
xlabel('whisker position');
fn = 'whisker_position_pop.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])


% onoff_heat = nan(2,numel(plot_units));
% [red_intersect] = ismember(sorted_idx,red_dots);
% [blue_intersect]=ismember(sorted_idx,blue_dots);
% [gray_intersect]=ismember(sorted_idx,gray_dots);
% 
% onoff_heat(1,red_intersect) = 1;
% onoff_heat(1,blue_intersect) = -1;
% onoff_heat(1,gray_intersect) = 0;
% 
% figure(46);clf
% pcolor(flipud([onoff_heat]'));
% colormap(redbluecmap)

% fn = 'whisker_position_pop_redblue.eps';
% export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
% fix_eps_fonts([saveDir, fn])

%% modulation depth of other features
selectedCells = 1:length(U);
% pole_whisk = whisking_location_quantification(U,selectedCells,'pole','off');
angle_whisk = whisking_location_quantification(U,selectedCells,'angle','off');
midpoint_whisk = whisking_location_quantification(U,selectedCells,'midpoint','off');
amp_whisk = whisking_location_quantification(U,selectedCells,'amplitude','off');
phase_whisk = whisking_location_quantification(U,selectedCells,'phase','off');

% w_mod_idx_pole = cellfun(@(x) x.calculations.mod_idx_relative,pole_whisk(selectedCells));
w_mod_idx_angle = cellfun(@(x) x.calculations.mod_idx_relative,angle_whisk(selectedCells));
w_mod_idx_phase = cellfun(@(x) x.calculations.mod_idx_relative,phase_whisk(selectedCells));
w_mod_idx_amp = cellfun(@(x) x.calculations.mod_idx_relative,amp_whisk(selectedCells));
w_mod_idx_midpoint = cellfun(@(x) x.calculations.mod_idx_relative,midpoint_whisk(selectedCells));

% w_peak_idx_pole = cellfun(@(x) x.calculations.tune_peak,pole_whisk(selectedCells)) * -1;
w_peak_idx_angle = cellfun(@(x) x.calculations.tune_peak,angle_whisk(selectedCells));
w_peak_idx_phase = cellfun(@(x) x.calculations.tune_peak,phase_whisk(selectedCells));
w_peak_idx_amp = cellfun(@(x) x.calculations.tune_peak,amp_whisk(selectedCells));
w_peak_idx_midpoint = cellfun(@(x) x.calculations.tune_peak,midpoint_whisk(selectedCells));

% HEAT GRAY WHISK SORTING
[~,idx] = sort(w_mod_idx_angle);
whisk_angle_tuned= cellfun(@(x) x.is_tuned==1,angle_whisk);
whisktunedidx = ismember(idx,find(whisk_angle_tuned));
nonwhisktunedidx = ismember(idx,find(~whisk_angle_tuned));
tune_map = [w_mod_idx_angle ; w_mod_idx_phase ; w_mod_idx_amp ; w_mod_idx_midpoint];
fr_map  =  firing_rate;

figure(29);clf
subplot(2,2,1)
imagesc(tune_map(:,idx(whisktunedidx)))
set(gca,'ytick',1:6,'yticklabel',{'angle','phase','amp','midpoint'})
caxis([0 1])
subplot(2,2,2)
imagesc(tune_map(:,idx(nonwhisktunedidx)))
set(gca,'ytick',1:6,'yticklabel',{'angle','phase','amp','midpoint'})
caxis([0 1])
subplot(2,2,3)
imagesc(fr_map(:,idx(whisktunedidx)))
set(gca,'ytick',1,'yticklabel',{'log firing rate'})
caxis([-2 2])
subplot(2,2,4)
imagesc(fr_map(:,idx(nonwhisktunedidx)))
set(gca,'ytick',1:2,'yticklabel',{'log firing rate'})
caxis([-2 2])
colormap(gray)

fn = 'gray_whisk_sort.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])

figure(3480);clf
subplot(1,2,1)
imagesc(abs(corr([tune_map(:,idx(whisktunedidx)) ; fr_map(:,idx(whisktunedidx))]')));
caxis([0 1])
set(gca,'xticklabel',{'angle','phase','amp','midpoint','firing rate'},'ytick',[]);xtickangle(45)
axis square 
subplot(1,2,2)
imagesc(abs(corr([tune_map(:,idx(nonwhisktunedidx)) ; fr_map(:,idx(nonwhisktunedidx))]')));
set(gca,'xticklabel',{'angle','phase','amp','midpoint','firing rate'},'ytick',[])
xtickangle(45)
colormap gray
caxis([0 1])
axis square 

fn = 'modulation_correlation.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])

% HSV PLOTTING WHISK
[~,sort_idx] = sort(w_peak_idx_angle);

peaks = {w_peak_idx_angle,w_peak_idx_phase,w_peak_idx_amp,w_peak_idx_midpoint};
mod_idx = {w_mod_idx_angle,w_mod_idx_phase,w_mod_idx_amp,w_mod_idx_midpoint};
final_image = nan(numel(peaks),numel(selectedCells),3);
for b = 1:numel(peaks)
    hues = normalize_var(peaks{b}(sort_idx),.7,1);
    hsv = [hues' mod_idx{b}(sort_idx)' ones(numel(hues),1)];
    rgb_values = hsv2rgb(hsv);
    final_image(b,:,1) = rgb_values(:,1);
    final_image(b,:,2) = rgb_values(:,2);
    final_image(b,:,3) = rgb_values(:,3);
end
%
figure(580);clf
imshow(final_image)
% saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig2\';
% fn = 'hsv_whisk.eps';
% export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
% fix_eps_fonts([saveDir, fn])
% 
% % HEAT GRAY WHISK TOUCH TUNING SORTING 
% [~,idx] = sort(w_mod_idx_angle);
% nontouchidx = ismember(idx,nontouchCells);
% tunedidx = ismember(idx,tunedCells);
% nontunedidx = ismember(idx,nontunedCells);
% 
% tune_map = [w_mod_idx_angle ; w_mod_idx_phase ; w_mod_idx_amp ; w_mod_idx_midpoint];
% fr_map  =  firing_rate;
% figure(19);clf
% subplot(2,3,1)
% imagesc(tune_map(:,idx(nontouchidx)))
% set(gca,'ytick',1:6,'yticklabel',{'angle','phase','amp','midpoint'})
% caxis([0 1])
% subplot(2,3,2)
% imagesc(tune_map(:,idx(tunedidx)))
% caxis([0 1])
% subplot(2,3,3)
% imagesc(tune_map(:,idx(nontunedidx)))
% caxis([0 1])
% colormap(gray)
% colorbar
% 
% subplot(2,3,4)
% imagesc(fr_map(:,idx(nontouchidx)))
% set(gca,'ytick',1,'yticklabel',{'log firing rate'})
% caxis([0 1])
% subplot(2,3,5)
% imagesc(fr_map(:,idx(tunedidx)))
% caxis([0 1])
% subplot(2,3,6)
% imagesc(fr_map(:,idx(nontunedidx)))
% caxis([-2 2])
% colormap(gray)
% colorbar
% 
% fn = 'gray_whisk_touch.eps';
% export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
% fix_eps_fonts([saveDir, fn])






