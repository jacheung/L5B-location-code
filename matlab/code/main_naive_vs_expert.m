%% Load whisking and neural time series struct
clear
load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\excitatory_all.mat') %L5b excitatory cells recorded by Jon and Phil
% load('C:\Users\jacheung\Dropbox\LocationCode\DataStructs\interneurons.mat') %L5b inhibitory cells

%% Top level parameters and definitions
hilbertVar = 'pole';
selectedCells = find(cellfun(@(x) strcmp(x.meta.touchProperties.responseType,'excited'),U));
pole_tuned = object_location_quantification(U,selectedCells,hilbertVar,'off'); %for old see object_location_v1.0
%% accuracy during recording sessions
expert = cellfun(@(x) strcmp(x.meta.layer,'BVL5b'),U);
naive = cellfun(@(x) ~strcmp(x.meta.layer,'BVL5b'),U);
recording_accuracy = cellfun(@(x) mean(x.meta.trialCorrect),U(expert));
naive_accuracy = cellfun(@(x) mean(x.meta.trialCorrect),U(naive));

figure(51);clf
scatter(ones(1,sum(naive)),naive_accuracy,'markeredgecolor',[.7 .7 .7]);
hold on; errorbar(1,mean(naive_accuracy),std(naive_accuracy),'ko'); 
scatter(ones(1,sum(expert)).*2,recording_accuracy,'markeredgecolor',[.7 .7 .7]);
hold on; errorbar(2,mean(recording_accuracy),std(recording_accuracy),'ko'); 
hold on; plot([0 3],[.5 .5],'--k')
set(gca,'ylim',[0 1],'ytick',0:.25:1','xtick',[],'xlim',[0 3])
% 
% saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig4\';
% fn = 'naive_expert_performance.eps';
% export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
% fix_eps_fonts([saveDir, fn])

rawPsychometricCurves(U(expert))

saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig2\';
fn = 'expert_psycho.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])
%% stimulus comparison of naive vs expert
% num touches
touch_nums = cellfun(@(x) sum(squeeze(x.S_ctk(9,:,:)==1) + squeeze(x.S_ctk(13,:,:)==1)),U,'uniformoutput',0);
whisk_mat = cellfun(@(x) squeeze(x.S_ctk(3,:,:)>5),U,'uniformoutput',0);

w_prop_expert = cellfun(@(x) mean(x(:)),whisk_mat(expert));
w_prop_naive = cellfun(@(x) mean(x(:)),whisk_mat(naive));

expert_mean = [cellfun(@(x) nanmean(x),touch_nums(expert))];
naive_mean = [cellfun(@(x) nanmean(x),touch_nums(naive))];
        
figure(48);clf
subplot(1,2,1)
histogram(expert_mean,0:1:20,'normalization','probability','facealpha',1,'facecolor','y')
hold on; histogram(naive_mean,0:1:20,'normalization','probability','facecolor',[.8 .8 .8])
set(gca,'xtick',0:5:20,'xlim',[0 20],'ylim',[0 .25],'ytick',0:.25:1)
xlabel('touch counts')
ylabel('proportion of units')
[~,p_counts] = ttest2(expert_mean,naive_mean);
title(num2str(p_counts));


figure(48);subplot(1,2,2)
histogram(w_prop_expert,0:.05:1,'normalization','probability','facealpha',1,'facecolor','y')
hold on; histogram(w_prop_naive,0:.05:1,'normalization','probability','facecolor',[.8 .8 .8])
set(gca,'xtick',0:.5:1,'xlim',[0 1],'ylim',[0 .5],'ytick',0:.25:1)
xlabel('proportion of time whisking')
ylabel('proportion of units')
[~,p_whisk] = ttest2(w_prop_expert,w_prop_naive);
title(num2str(p_whisk));

saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig2\';
fn = 'naive_expert_stimulus.eps';
export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
fix_eps_fonts([saveDir, fn])


%% Comparison of naive vs expert proportion of touch/OL cells
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
subplot(1,2,1)
hold on; bar(1:2,[naive_proportion_touch expert_proportion_touch],'facecolor',[.8 .8 .8])
hold on; bar(1:2,[naive_proportion_OL expert_proportion_OL],'facecolor','k')
ylabel('proportion of units')
set(gca,'xtick',1:2,'xticklabel',{['naive n=' num2str(sum(naive))],['expert n=' num2str(sum(expert))]},'ytick',0:.25:1,'ylim',[0 1])
legend('touch','location tuned','location','northwest')

subplot(1,2,2)
hold on; bar(1:2,[num_naive_OL./num_naive_touch, num_expert_OL./num_expert_touch],'facecolor','k')

ylabel('proportion of units')
set(gca,'xtick',1:2,'xticklabel',{['naive n=' num2str(num_naive_OL)],['expert n=' num2str(num_expert_OL)]},'ytick',0:.25:1,'ylim',[0 1])
title('proportion of touch that is OL')

    figure(99)
    saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig2\';
    fn = 'unit_distribution_bar.eps';
    export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
    fix_eps_fonts([saveDir, fn])

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
    unsorted_heat = norm_new(cell2mat(touch_heat')');
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
colormap turbo


    figure(80);
    saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig2\';
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
        norm_y = norm_new(sr.values(:,2));
        
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
    
    set(gca,'xlim',[-2 2],'xdir','reverse','ytick',0:.5:1,'xtick',-2:1:2)
    xlabel('distance from peak (mm)')
    axis square
    
end
       
   
        figure(30);
    saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig3\';
    fn = 'naive_expert_modulation_width.eps';
    export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
    fix_eps_fonts([saveDir, fn])
    



