%% Load/build touch structures 
clear
data_directory = 'C:\Users\jacheung\Dropbox\LocationCode\DataStructs\';
load([data_directory 'Raw\excitatory_all.mat']);
feature_list = {'pole'};
touch_struct = CompileTStruct(data_directory, feature_list);
%% Raster + PSTH one cell (A) 
trial_number = 29; %example trial in publication
plot_trial_raster(U, trial_number);
%% firing rate X depth of recording (B)
plot_cell_depth(U, touch_struct)

%% touch psth by quartiles of far, close and near (C)
units_to_plot = [13 ,5, 8];
plot_example_PSTH(U, touch_struct, units_to_plot);

%% heatmap for object location tuned touch units (D)
variable = 'pole';
tuned_structs = pole_tuned(cellfun(@(x) x.is_tuned==1,pole_tuned));
touch_heat = cell(1,numel(tuned_structs));

for g = 1:numel(tuned_structs)
    curr_t = tuned_structs{g}.stim_response.values;
    curr_t = curr_t(~any(isnan(curr_t),2),:);%clean nan rows
    [~,u_idx] = unique(curr_t(:,1)); %catch non-unique x-values
    curr_t = curr_t(u_idx,:);
    if strcmp(variable,'pole')
        touch_x = -1:.1:1;
    elseif strcmp(variable,'angle')
        touch_x = linspace(-30,60,21);
    end
    
    touch_heat{g} = interp1(curr_t(:,1),curr_t(:,2),touch_x);
end
unsorted_heat = norm_new(cell2mat(touch_heat')');
[~,t_max_idx] = max(unsorted_heat,[],1);
[~,t_idx] = sort(t_max_idx);
data = unsorted_heat(:,t_idx)';
%set unsampled heatmap to nan;
[nr,nc] = size(data);
figure(15);clf
pcolor([data nan(nr,1); nan(1,nc+1)]);
caxis([0 1])
shading flat;
colormap turbo
set(gca,'xdir','reverse','xtick',1:10:length(touch_x),'xlim',[1 length(touch_x)],'xticklabel',-1:1:1,'ydir','reverse')
title(['location tuned heat map (n=' num2str(size(data,1)) ')'])
colorbar

fn = 'heat_map_tuning.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])

%% modulation width (E)
% pole_tuned = object_location_quantification(U,selectedCells,'pole','off'); %for old see object_location_v1.0
tuned_units = find(cellfun(@(x) x.is_tuned==1,pole_tuned));

peak_response = cellfun(@(x) x.calculations.tune_peak,pole_tuned(tuned_units));
interp_norm_y = cell(1,length(peak_response));
for g = 1:length(peak_response)
    sr = pole_tuned{tuned_units(g)}.stim_response;
    centered_x = sr.values(:,1) - peak_response(g) ;
    norm_y = norm_new(sr.values(:,2));
    
    interp_centered_x = -2:.1:2;
    raw_x = round(centered_x,2);
    [~,idx] = unique(raw_x);
    interp_norm_y{g} = interp1(raw_x(idx),norm_y(idx),interp_centered_x);
    figure(30);
    hold on; plot(interp_centered_x,interp_norm_y{g},'color',[.8 .8 .8])
end

pop_mean = nanmean(cell2mat(interp_norm_y'));
pop_sem = nanstd(cell2mat(interp_norm_y')) ./ sqrt(sum(~isnan(cell2mat(interp_norm_y'))));
hold on; shadedErrorBar(interp_centered_x,pop_mean,pop_sem,'k')

set(gca,'xlim',[-2 2],'xdir','reverse','ytick',0:.5:1,'xtick',-2:1:2)
xlabel('distance from peak (mm)')
axis square

figure(30);
fn = 'all_modulation_width.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])

%% Follicle at first touch vs later (SF)
pxpermm = 33;
clear fdist_map
for i = 1:numel(U)
    for b = 1:U{i}.k
        touchOn = [find(U{i}.S_ctk(9,:,b)==1) find(U{i}.S_ctk(12,:,b)==1)];
        if numel(touchOn)>1
            fx = U{i}.whisker.follicleX{b};
            fy = U{i}.whisker.follicleY{b};
            touchOn = touchOn(touchOn < numel(fx));
            fxy = [fx(touchOn)' fy(touchOn)'];
            fdist_map_x{i}{b} = fxy(2:end,1) - fxy(1,1);
            fdist_map_y{i}{b} = fxy(2:end,2) - fxy(1,2);
%             fdist = pdist2(fxy,fxy);
%             fdist_map{i}{b} = diag(fdist,1)';
        end
    end
end

num_touches_per_trial = cellfun(@(x) cellfun(@(y) numel(y),x), fdist_map_x,'uniformoutput',0); 
touch_num_cdf = cellfun(@(x) cumsum(histc(x,0:50)./numel(x)),num_touches_per_trial,'uniformoutput',0);
indiv = cell2mat(touch_num_cdf');
pop = mean(indiv);
pop_num_sd = std(indiv);
pop__num_sem = std(indiv)./sqrt(numel(U)); 

figure(80);clf
subplot(1,3,1)
plot(0:50,indiv,'color',[.8 .8 .8])
hold on;shadedErrorBar(0:50,pop,pop_num_sd,'r');
set(gca,'ylim',[0 1],'ytick',0:.25:1)
xlabel('number of touches per trial')
ylabel('proportion of all trials')

all_mean_fd = cellfun(@(x) nanmean(cell2nanmat(x),2), fdist_map_x,'uniformoutput',0);
pop_mean = nanmean(cell2nanmat(all_mean_fd),2) * -1;
pop_sd = nanstd(cell2nanmat(all_mean_fd),[],2);
pop_sem = nanstd(cell2nanmat(all_mean_fd),[],2) ./ sqrt(sum(~isnan(cell2nanmat(all_mean_fd)),2));
subplot(1,3,2)
errorbar(2:numel(pop_mean)+1,pop_mean./pxpermm,pop_sd./pxpermm,'-ko')
set(gca,'xlim',[0 20],'ylim',[-.5 .5],'ytick',-.5:.25:.5)
ylabel('anterior --- posterior')
xlabel('touch number')

all_mean_fd = cellfun(@(x) nanmean(cell2nanmat(x),2), fdist_map_y,'uniformoutput',0);
pop_mean = nanmean(cell2nanmat(all_mean_fd),2);
pop_sd = nanstd(cell2nanmat(all_mean_fd),[],2);
pop_sem = nanstd(cell2nanmat(all_mean_fd),[],2) ./ sqrt(sum(~isnan(cell2nanmat(all_mean_fd)),2));
subplot(1,3,3)
errorbar(2:numel(pop_mean)+1,pop_mean./pxpermm,pop_sd./pxpermm,'-ko')
set(gca,'xlim',[0 20],'ylim',[-.5 .5],'ytick',-.5:.25:.5)
ylabel('medial --- lateral')
xlabel('touch number')

saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig3\';
fn = 'scatter_depth_firingrate.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])







