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
%% heatmap
population_heatmap_builder(tStruct,wStruct,hilbertVar)

%% proportion of neurons that intersect bar chart
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
saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig6\';
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

%% intersection of whisking and touch
tUnits = cellfun(@(x) x.is_tuned==1,tStruct);
wUnits = cellfun(@(x) x.is_tuned==1,wStruct);

[~,touch_ix_idx] = intersect(find(tUnits),find(wUnits));
[~,whisk_ix_idx] = intersect(find(wUnits),find(tUnits));

touch_nonIX_idx = setdiff(1:sum(tUnits),touch_ix_idx);
whisk_nonIX_idx = setdiff(1:sum(wUnits),whisk_ix_idx);

whiskTuned = find(wUnits);
touchTuned = find(tUnits);
touch_whisk_tuned = intersect(find(tUnits),find(wUnits));

touch_OL = logical(ones(1,length(U)));
rc = numSubplots(sum(touch_OL));

sel_tstructs = tStruct(touch_OL);
sel_wstructs = wStruct(touch_OL);

figure(100);clf
figure(101);clf
whisk_touch_pair = cell(1,sum(touch_OL));
touch_diff_pair = cell(1,sum(touch_OL));
for g = 1:sum(touch_OL)
    if isfield(sel_wstructs{g},'stim_response') && isfield(sel_tstructs{g},'stim_response')
        curr_w = sel_wstructs{g}.stim_response.values;
        curr_t = sel_tstructs{g}.stim_response.values;
        
        curr_w = curr_w(~any(isnan(curr_w),2),:); %clean nan rows
        curr_t = curr_t(~any(isnan(curr_t),2),:);
        
        if strcmp(hilbertVar,'pole')
            whisk_x = round(round(min(curr_w(:,1)),1):.1:round(max(curr_w(:,1)),1),1);
            touch_x = round(round(min(curr_t(:,1)),1):.1:round(max(curr_t(:,1)),1),1);
        elseif strcmp(hilbertVar,'phase')
            whisk_x = linspace(-pi,pi,21);
            touch_x = linspace(-pi,pi,21);
        else
            whisk_x = round(round(min(curr_w(:,1))):1:round(max(curr_w(:,1))));
            touch_x = round(round(min(curr_t(:,1))):1:round(max(curr_t(:,1))));
        end
        whisk_response = interp1(curr_w(:,1),curr_w(:,2),whisk_x);
        whisk_std = interp1(curr_w(:,1),curr_w(:,3),whisk_x);
        whisk_CI = interp1(curr_w(:,1),curr_w(:,4),whisk_x);
        
        touch_response = interp1(curr_t(:,1),curr_t(:,2),touch_x);
        touch_std = interp1(curr_t(:,1),curr_t(:,3),touch_x);
        touch_CI =  interp1(curr_t(:,1),curr_t(:,4),touch_x);
        
        [~,~,whisk_idx] = intersect(touch_x,whisk_x);
        [overlap_x,~,touch_idx] = intersect(whisk_x,touch_x);
        
        %raw responses within touch ranges
        %     figure(99);subplot(rc(1),rc(2),g)
        %     shadedErrorBar(overlap_x,whisk_response(whisk_idx),whisk_CI(whisk_idx),'c')
        %     hold on; shadedErrorBar(overlap_x,touch_response(touch_idx),touch_CI(touch_idx),'b')
        %     if strcmp(hilbertVar,'pole')
        %         set(gca,'xlim',[-1 1],'xdir','reverse')
        %     elseif strcmp(hilbertVar,'phase')
        %         set(gca,'xlim',[-pi pi],'xtick',-pi:pi:pi,'xticklabel',{'\pi','0','\pi'})
        %     end
        
        %raw responses
        figure(100);subplot(rc(1),rc(2),g)
        shadedErrorBar(whisk_x,whisk_response,whisk_CI,'c')
        hold on; shadedErrorBar(touch_x,touch_response,touch_CI,'b')
        if strcmp(hilbertVar,'pole')
            set(gca,'xlim',[-1 2],'xdir','reverse')
            axis square
        elseif strcmp(hilbertVar,'phase')
            set(gca,'xlim',[-pi pi],'xtick',-pi:pi:pi,'xticklabel',{'\pi','0','\pi'})
        end
        
        if any(g == whiskTuned(whisk_nonIX_idx))
            title('whisk only')
        elseif any(g == touchTuned(touch_nonIX_idx))
            title('touch only')
        elseif any(g ==touch_whisk_tuned)
            title('T+W')
        end
        
        %normalized responses MAXMIN
%         m = min(whisk_response);
%         range = max(whisk_response) - m;
%         norm_whisk= (whisk_response - m) ./ range;
%         norm_whisk_CI = (whisk_CI - m) ./ range;
%         
%         m = min(touch_response);
%         range = max(touch_response) - m;
%         norm_touch= (touch_response - m) ./ range;
%         norm_touch_CI = (touch_CI - m) ./ range;
%         
        %zscoring
        mu = nanmean(whisk_response);
        sigma = nanstd(whisk_response);
        norm_whisk= (whisk_response - mu) ./ sigma;
        norm_whisk_CI = (whisk_CI - mu) ./ sigma;
        
        mu = nanmean(touch_response);
        sigma = nanstd(touch_response);
        norm_touch= (touch_response - mu) ./ sigma;
        norm_touch_CI = (touch_CI - mu) ./ sigma;
        

    
        figure(101);subplot(rc(1),rc(2),g)
            shadedErrorBar(whisk_x,norm_whisk,norm_whisk_CI,'c')
            hold on; shadedErrorBar(touch_x,norm_touch,norm_touch_CI,'b')
%         plot(whisk_x,norm_whisk,'c')
%         hold on; plot(touch_x,norm_touch,'b')
        if strcmp(hilbertVar,'pole')
            set(gca,'xlim',[-1 2],'xdir','reverse','ylim',[-3 3])
            axis square
        elseif strcmp(hilbertVar,'phase')
            set(gca,'xlim',[-pi pi],'xtick',-pi:pi:pi,'xticklabel',{'\pi','0','\pi'})
        end
        
        if any(g == whiskTuned(whisk_nonIX_idx))
            title('whisk only')
        elseif any(g == touchTuned(touch_nonIX_idx))
            title('touch only')
        elseif any(g ==touch_whisk_tuned)
            title('T+W')
        end
        
    end
    
    
    
end

%     figure(100);
%     saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig6\';
%     fn = 'whisk_touch_tuning_curves.eps';
%     export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
%     fix_eps_fonts([saveDir, fn])

%     figure(101);
%     saveDir = 'C:\Users\jacheung\Dropbox\LocationCode\Figures\Parts\Fig5\';
%     fn = 'whisk_touch_tuning_curves_normalized.eps';
%     export_fig([saveDir, fn], '-depsc', '-painters', '-r1200', '-transparent')
%     fix_eps_fonts([saveDir, fn])

