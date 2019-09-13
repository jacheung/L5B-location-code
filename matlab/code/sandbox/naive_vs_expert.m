%% Load whisking and neural time series struct
clear
load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\excitatory_all.mat') %L5b excitatory cells recorded by Jon and Phil
% load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\interneurons.mat') %L5b inhibitory cells

%% Top level parameters and definitions
hilbertVar = 'pole'
selectedCells = find(cellfun(@(x) strcmp(x.meta.touchProperties.responseType,'excited'),U));
pole_tuned = object_location_quantification(U,selectedCells,hilbertVar,'off'); %for old see object_location_v1.0

%% Comparison of naive vs expert
naive = cellfun(@(x) ~strcmp(x.meta.layer,'BVL5b'),U);
expert = cellfun(@(x) strcmp(x.meta.layer,'BVL5b'),U);

num_naive_touch = sum(cellfun(@(x) strcmp(x.meta.touchProperties.responseType,'excited'),U(naive)));
num_expert_touch = sum(cellfun(@(x) strcmp(x.meta.touchProperties.responseType,'excited'),U(expert)));
num_naive_OL = sum((cellfun(@(x) x.is_tuned==1,pole_tuned(naive)))) ; 
num_expert_OL = sum((cellfun(@(x) x.is_tuned==1,pole_tuned(expert)))) ; 

naive_proportion_touch = num_naive_touch ./ sum(naive);
expert_proportion_touch = num_expert_touch ./ sum(expert);
naive_proportion_OL  = num_naive_OL ./ sum(naive);
expert_proportion_OL  = num_expert_OL ./ sum(expert);

figure(99);clf
hold on; bar(1:2,[naive_proportion_touch expert_proportion_touch],'facecolor','k')
hold on; bar(1:2,[naive_proportion_OL expert_proportion_OL],'facecolor',[.8 .8 .8])
ylabel('proportion of units')
set(gca,'xtick',1:2,'xticklabel',{['naive n=' num2str(sum(naive))],['expert n=' num2str(sum(expert))]},'ytick',0:.25:1,'ylim',[0 .75])
legend('touch','location tuned','location','northwest')
%     figure(99)
%     saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig3\';
%     fn = 'unit_distribution_bar.eps';
%     export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
%     fix_eps_fonts([saveDir, fn])

%% Heatmap of tuning
tuned_units = cellfun(@(x) x.is_tuned==1,pole_tuned);
naive = cellfun(@(x) ~strcmp(x.meta.layer,'BVL5b'),U);
expert = cellfun(@(x) strcmp(x.meta.layer,'BVL5b'),U);

naive_tuned_units = intersect(find(tuned_units),find(naive));
expert_tuned_units = intersect(find(tuned_units),find(expert));

tuned_units = {naive_tuned_units,expert_tuned_units};

figure(80);clf
for b=1:length(tuned_units)
    sel_tstructs = pole_tuned(tuned_units{b});
    
    touch_heat = cell(1,length(sel_tstructs));
    %building heatmap for object location tuned touch units
    for g = 1:numel(sel_tstructs)
        curr_t = sel_tstructs{g}.stim_response.values;
        
        %clean nan rows
        curr_t = curr_t(~any(isnan(curr_t),2),:);
        
        if strcmp(hilbertVar,'pole')
            touch_x = -1:.1:1;
        elseif strcmp(hilbertVar,'phase')
            touch_x = linspace(-pi,pi,21);
        else
            touch_x = round(round(min(curr_t(:,1))):1:round(max(curr_t(:,1))));
        end
        touch_heat{g} = interp1(curr_t(:,1),curr_t(:,2),touch_x);
    end
    
    %plotting all touch object location tuned units
    subplot(1,2,b)
    unsorted_heat = normalize_var(cell2mat(touch_heat')',0,1);
    [~,t_max_idx] = max(unsorted_heat,[],1);
    [~,t_idx] = sort(t_max_idx);
    data = unsorted_heat(:,t_idx)';
    %set unsampled heatmap to nan;
    [nr,nc] = size(data);
    pcolor([data nan(nr,1); nan(1,nc+1)]);
    shading flat;
    set(gca,'xdir','reverse','xtick',1:10:length(touch_x),'xlim',[1 length(touch_x)],'xticklabel',-1:1:1,'ydir','reverse')
    axis square
end


    figure(80);
    saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig3\';
    fn = 'naive_expert_heat.eps';
    export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
    fix_eps_fonts([saveDir, fn])
    
%% shape of tuning
tuned_units = cellfun(@(x) x.is_tuned==1,pole_tuned);
naive = cellfun(@(x) ~strcmp(x.meta.layer,'BVL5b'),U);
expert = cellfun(@(x) strcmp(x.meta.layer,'BVL5b'),U);

naive_tuned_units = intersect(find(tuned_units),find(naive));
expert_tuned_units = intersect(find(tuned_units),find(expert));

tuned_units_list = {naive_tuned_units,expert_tuned_units};

figure(30);clf
for b=1:length(tuned_units_list)
    peak_response = cellfun(@(x) x.calculations.tune_peak,pole_tuned(tuned_units_list{b}));
    for g = 1:length(peak_response)
        sr = pole_tuned{tuned_units_list{b}(g)}.stim_response;
        centered_x = sr.values(:,1) - peak_response(g) ;
        norm_y = normalize_var(sr.values(:,2),0,1);
        
        interp_centered_x = -2:.1:2;
        raw_x = round(centered_x,2);
        [~,idx] = unique(raw_x);
        interp_norm_y{g} = interp1(raw_x(idx),norm_y(idx),interp_centered_x);
        figure(30);subplot(1,2,b)
        hold on; plot(interp_centered_x,interp_norm_y{g},'color',[.8 .8 .8])
    end
    
    pop_mean = nanmean(cell2mat(interp_norm_y'));
    pop_sem = nanstd(cell2mat(interp_norm_y')) ./ sqrt(sum(~isnan(cell2mat(interp_norm_y'))));
    hold on; shadedErrorBar(interp_centered_x,pop_mean,pop_sem,'k')
    
    set(gca,'xlim',[-1 1],'xdir','reverse','ytick',0:.5:1,'xtick',-2:1:2)
    xlabel('distance from peak (mm)')
    axis square
    
end
       
   
        figure(30);
    saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig3\';
    fn = 'naive_expert_modulation_width.eps';
    export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
    fix_eps_fonts([saveDir, fn])
    



