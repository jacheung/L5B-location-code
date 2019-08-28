% selectedCells = find(cellfun(@(x) isfield(x.meta.touchProperties,'responseWindow'),U)~=0);

%% 
hilbertVar = 'pole'


selectedCells = find(cellfun(@(x) strcmp(x.meta.touchProperties.responseType,'excited'),U));
tStruct = object_location_quantification(U,selectedCells,hilbertVar,'off');

wStruct = whisking_location_quantification(U,selectedCells,hilbertVar,'off');
%% intersection of whisking and touch 

touch_OL = cellfun(@(x) x.is_tuned==1,tStruct);
rc = numSubplots(sum(touch_OL)); 

sel_tstructs = tStruct(touch_OL);
sel_wstructs = wStruct(touch_OL);

figure(8);clf
figure(99);clf
whisk_touch_pair = cell(1,sum(touch_OL));
touch_diff_pair = cell(1,sum(touch_OL)); 
for g = 1:sum(touch_OL)
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
    
    figure(99);subplot(rc(1),rc(2),g)
    shadedErrorBar(overlap_x,whisk_response(whisk_idx),whisk_CI(whisk_idx),'c')
    hold on; shadedErrorBar(overlap_x,touch_response(touch_idx),touch_CI(touch_idx),'b')
    if strcmp(hilbertVar,'pole')
        set(gca,'xlim',[-1 1],'xdir','reverse')
    elseif strcmp(hilbertVar,'phase')
        set(gca,'xlim',[-pi pi],'xtick',-pi:pi:pi,'xticklabel',{'\pi','0','\pi'})
    else
        
    end
    
    
    figure(100);subplot(rc(1),rc(2),g)
    scatter(whisk_response(whisk_idx),touch_response(touch_idx))
    
    
%     plot(overlap_x,whisk_response(whisk_idx),'c')
%     hold on; plot(overlap_x,touch_response(touch_idx),'b')
    all_responses = [whisk_response(whisk_idx) touch_response(touch_idx)];
%     set(gca,'xlim',[min(all_responses) max(all_responses)],'ylim',[min(all_responses) max(all_responses)])
%     hold on; plot([min(all_responses) max(all_responses)],[min(all_responses) max(all_responses)],'-.k')
    
    normed_responses = normalize_var(all_responses,0,1);
    whisk_touch_pair{g} = reshape(normed_responses,numel(normed_responses)./2,2); 
    
    %calculation of difference to quantify the effect touch has on whisking
    response_difference = touch_response(touch_idx) - whisk_response(whisk_idx);
    figure(8);subplot(rc(1),rc(2),g)
    scatter(touch_response(touch_idx),response_difference,'b')
    set(gca,'xlim',[min([all_responses response_difference]) max([all_responses response_difference])],'ylim',[min([all_responses response_difference]) max([all_responses response_difference])])
    axis square
    hold on;plot([min([all_responses response_difference]) max([all_responses response_difference])],[min([all_responses response_difference]) max([all_responses response_difference])],'-.k')
    normed_responses_tdpair = normalize_var([response_difference touch_response(touch_idx)],0,1);
    touch_diff_pair{g} =  reshape(normed_responses_tdpair,numel(normed_responses_tdpair)./2,2); 

    
end   

all_values = cell2mat(whisk_touch_pair');
figure(12);clf
scatter(all_values(:,1),all_values(:,2),'k')
hold on; plot([0 1],[0 1],'-.k')
set(gca,'xlim',[0 1],'ylim',[0 1])
ylabel('normalized touch responses')
xlabel('normalized whisk responses')
axis square

pair_mean = cell2mat(cellfun(@nanmean,whisk_touch_pair','uniformoutput',0));
pair_sem = cell2mat(cellfun(@(x) nanstd(x)./sqrt(length(x)),whisk_touch_pair,'uniformoutput',0)');
figure(15);clf
subplot(1,2,1)
for g = 1:size(pair_mean,1)
    hold on;errorbar(pair_mean(g,1),pair_mean(g,2),pair_sem(g,2),pair_sem(g,2),pair_sem(g,1),pair_sem(g,1),'ko')
end
hold on; plot([0 1],[0 1],'--k')
set(gca,'xlim',[0 1],'ylim',[0 1])
axis square
ylabel('normalized touch responses')
xlabel('normalized whisk responses')

% all_values_tdpair = cell2mat(touch_diff_pair');
% figure(13);clf
% scatter(all_values_tdpair(:,2),all_values_tdpair(:,1),'k')
% hold on; plot([0 1],[0 1],'-.k')
% set(gca,'xlim',[0 1],'ylim',[0 1])
% ylabel('normalized touch responses')
% xlabel('normalized difference (touch-whisking)')
% axis square
% lm = fitlm(all_values_tdpair(:,1),all_values_tdpair(:,2))

pair_mean = cell2mat(cellfun(@nanmean,touch_diff_pair','uniformoutput',0));
pair_sem = cell2mat(cellfun(@(x) nanstd(x)./sqrt(length(x)),touch_diff_pair,'uniformoutput',0)')
subplot(1,2,2)
for g = 1:size(pair_mean,1)
    hold on;errorbar(pair_mean(g,2),pair_mean(g,1),pair_sem(g,1),pair_sem(g,1),pair_sem(g,2),pair_sem(g,2),'ko')
end
hold on; plot([0 1],[0 1],'--k')
set(gca,'xlim',[0 1],'ylim',[0 1])
axis square
xlabel('normalized touch responses')
ylabel('normalized difference (touch-whisking)')
suptitle(hilbertVar)



    
