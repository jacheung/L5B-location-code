saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig4\';
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

% figure(480);clf
% scatter(fr_quiet(red_dots),fr_whisk(red_dots),'filled','r')
% hold on; scatter(fr_quiet(blue_dots),fr_whisk(blue_dots),'filled','b')
% hold on; scatter(fr_quiet(gray_dots),fr_whisk(gray_dots),'filled','markerfacecolor',[.8 .8 .8])
% axis square
% hold on; plot([0 max([fr_quiet fr_whisk])],[0 max([fr_quiet fr_whisk])],'--k')
% set(gca,'xlim',[0 max([fr_quiet fr_whisk])],'ylim',[0 max([fr_quiet fr_whisk])])
% xlabel('quiet FR');ylabel('whisking FR')
% title(['red=' num2str(numel(red_dots)) ' blue=' num2str(numel(blue_dots)) ' gray=' num2str(numel(gray_dots))])
% 


%% tuned units

selectedCells = 1:length(U);
angle_whisk = whisking_location_quantification(U,selectedCells,'angle','off');

tuned_units = cellfun(@(x) x.is_tuned,angle_whisk)==1;



%% pop map 
built = find(cellfun(@(x) isfield(x,'stim_response'),angle_whisk)); 
bars_built = cellfun(@(x) isfield(x.stim_response,'bars_fit'),angle_whisk(built));
full_build = built(bars_built);

% plot_units = intersect(full_build);
plot_units = full_build

bars_io = cellfun(@(x) [x.stim_response.bars_stim' x.stim_response.bars_fit.mean ],angle_whisk(plot_units),'uniformoutput',0);

whisk_sampled = round([min(cellfun(@(x) min(x(:,1)),bars_io)):max(cellfun(@(x) max(x(:,1)),bars_io))]);

interped_resp = cellfun(@(x) interp1(x(:,1),x(:,2),whisk_sampled),bars_io,'uniformoutput',0);
norm_resp = cellfun(@(x) normalize_var(x,0,1),interped_resp,'uniformoutput',0);

[~,maxidx] = cellfun(@(x) max(x),norm_resp);
[~,sorted_idx] = sort(maxidx); 

sort_whisk_tuning = cell2mat(norm_resp(sorted_idx)'); 

figure(48);clf
for g = 1:numel(sorted_idx)
    hold on; plot(sort_whisk_tuning(g,:)+g,'k')
end


figure(45);clf
pcolor(flipud(cell2mat(norm_resp(sorted_idx)'))); %flipud because pcolor flips compared to imagesc
set(gca,'xtick',1:10:numel(whisk_sampled),'xticklabel',whisk_sampled(1:10:end))
xlabel('whisker position');
fn = 'whisker_position_pop.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])


onoff_heat = nan(2,numel(plot_units));
[red_intersect] = ismember(sorted_idx,red_dots);
[blue_intersect]=ismember(sorted_idx,blue_dots);
[gray_intersect]=ismember(sorted_idx,gray_dots);

onoff_heat(1,red_intersect) = 1;
onoff_heat(1,blue_intersect) = -1;
onoff_heat(1,gray_intersect) = 0;

tuned_heat = ismember(sorted_idx,find(tuned_units));

figure(46);clf
pcolor(flipud([tuned_heat;onoff_heat]'));
colormap(redbluecmap)


fn = 'whisker_position_pop_redblue.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])


%% Proportion of units tuned 



    
        saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig4\';
        fn = 'whisker_position_indiv.eps';
        export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
        fix_eps_fonts([saveDir, fn])