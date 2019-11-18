function population_heatmap_builder(tStruct,wStruct,hilbertVar)
% this function plots heatmaps of whisking and touch tuning 
%
% requires products of object_location_quantification and
% whisking_location_quantification first before running. 
%
% uses touch tuning struct and whisk tuning struct for POLE position tuning
% to build population heatmaps 

%% Plot heatmap of tuning across all units
tuned_units = cellfun(@(x) x.is_tuned==1,tStruct);
whisk_units = cellfun(@(x) x.is_tuned==1,wStruct);

sel_tstructs = tStruct(tuned_units);
sel_wstructs = wStruct(whisk_units); 

touch_heat = cell(1,sum(tuned_units));
whisk_heat = cell(1,sum(whisk_units)); 

touch_ix_tmp = intersect(find(tuned_units),find(whisk_units)); 
whisk_ix_tmp = intersect(find(whisk_units),find(tuned_units)); 

[~,touch_ix_idx] = intersect(find(tuned_units),touch_ix_tmp);
[~,whisk_ix_idx] = intersect(find(whisk_units),whisk_ix_tmp);

%building heatmap for object location tuned touch units
for g = 1:numel(sel_tstructs)
    curr_t = sel_tstructs{g}.stim_response.values;
    
    %clean nan rows
    curr_t = curr_t(~any(isnan(curr_t),2),:); 
    
    if strcmp(hilbertVar,'pole')
        touch_x = -1:.1:1;
    elseif strcmp(hilbertVar,'phase')
        touch_x = linspace(-pi,pi,21);
    elseif strcmp(hilbertVar,'angle')
        touch_x = -30:80;
    else
        touch_x = round(round(min(curr_t(:,1))):1:round(max(curr_t(:,1))));
    end
    touch_heat{g} = interp1(curr_t(:,1),curr_t(:,2),touch_x);
    
end

%building heatmap for object location tuned whisk units
for d = 1:numel(sel_wstructs)
    curr_w = sel_wstructs{d}.stim_response.values;
    
    %clean nan rows
    curr_w = curr_w(~any(isnan(curr_w),2),:);
    
    if strcmp(hilbertVar,'pole')
        whisk_x = -2:.1:2;
    elseif strcmp(hilbertVar,'phase')
        whisk_x = linspace(-pi,pi,21);
    elseif strcmp(hilbertVar,'angle')
        whisk_x = -30:80;
    else
        whisk_x = round(round(min(curr_t(:,1))):1:round(max(curr_t(:,1))));
    end
    whisk_heat{d} = interp1(curr_w(:,1),curr_w(:,2),whisk_x);
end
%%
figure(50);clf

%plotting all touch object location tuned units
subplot(3,2,1)
% unsorted_heat = normalize_var(cell2mat(touch_heat')',0,1);
unsorted_heat = norm_new(cell2mat(touch_heat')');
[~,t_max_idx] = max(unsorted_heat,[],1);
[~,t_idx] = sort(t_max_idx);
data = unsorted_heat(:,t_idx)'; 
%set unsampled heatmap to nan; 
[nr,nc] = size(data);
pcolor([data nan(nr,1); nan(1,nc+1)]);
caxis([0 1])
shading flat;
colormap turbo
if strcmp(hilbertVar,'pole')
    set(gca,'xdir','reverse','xtick',1:10:length(touch_x),'xlim',[1 length(touch_x)],'xticklabel',-1:1:1,'ydir','reverse')
elseif strcmp(hilbertVar,'phase')
    set(gca,'xtick',1:10:length(touch_x),'xlim',[1 length(touch_x)],'xticklabel',-1:1:1,'ydir','reverse')
elseif strcmp(hilbertVar,'angle')
    set(gca,'xtick',10:20:length(whisk_x),'xticklabel',-20:20:80,'xlim',[0 length(whisk_x)])
end
title('all touch sorted')

%plotting all whisk object location tuned units
subplot(3,2,2)
% unsorted__whisk_heat = normalize_var(cell2mat(whisk_heat')',0,1);
unsorted__whisk_heat = norm_new(cell2mat(whisk_heat')');
[~,w_max_idx] = max(unsorted__whisk_heat,[],1);
[~,w_idx] = sort(w_max_idx);
data = unsorted__whisk_heat(:,w_idx)'; 
%set unsampled heatmap to nan; 
[nr,nc] = size(data);
pcolor([data nan(nr,1); nan(1,nc+1)]);
caxis([0 1])
shading flat;
colormap turbo
if strcmp(hilbertVar,'pole')
    set(gca,'xdir','reverse','xtick',1:10:length(touch_x),'xlim',[1 length(touch_x)],'xticklabel',-1:1:1,'ydir','reverse')
elseif strcmp(hilbertVar,'phase')
    set(gca,'xtick',1:10:length(touch_x),'xlim',[1 length(touch_x)],'xticklabel',-1:1:1,'ydir','reverse')
elseif strcmp(hilbertVar,'angle')
    set(gca,'xtick',10:20:length(whisk_x),'xticklabel',-20:20:80,'xlim',[0 length(whisk_x)])
end
title('all whisk sorted')

%plotting all whisk+touch object location tuned units touch responses sorted by touch peak 
subplot(3,2,3)
% unsorted_heat_touch = normalize_var(cell2mat(touch_heat(touch_ix_idx)')',0,1);
unsorted_heat_touch = norm_new(cell2mat(touch_heat(touch_ix_idx)')');
[~,t_max_idx] = max(unsorted_heat_touch,[],1);
[~,t_idx_ix] = sort(t_max_idx);
data = unsorted_heat_touch(:,t_idx_ix)'; 
%set unsampled heatmap to nan; 
[nr,nc] = size(data);
pcolor([data nan(nr,1); nan(1,nc+1)]);
caxis([0 1])
shading flat;
colormap turbo
if strcmp(hilbertVar,'pole')
    set(gca,'xdir','reverse','xtick',1:10:length(touch_x),'xlim',[1 length(touch_x)],'xticklabel',-1:1:1,'ydir','reverse')
elseif strcmp(hilbertVar,'phase')
    set(gca,'xtick',1:10:length(touch_x),'xlim',[1 length(touch_x)],'xticklabel',-1:1:1,'ydir','reverse')
elseif strcmp(hilbertVar,'angle')
    set(gca,'xtick',10:20:length(whisk_x),'xticklabel',-20:20:80,'xlim',[0 length(whisk_x)])
end
title('ix touch sorted')

%plotting all whisk+touch object location tuned units whisk responses sorted by whisk peak 
subplot(3,2,4)
% unsorted_whisk_heat = normalize_var(cell2mat(whisk_heat(whisk_ix_idx)')',0,1);
unsorted_whisk_heat = norm_new(cell2mat(whisk_heat(whisk_ix_idx)')');
[~,w_max_idx] = max(unsorted_whisk_heat,[],1);
[~,w_idx_ix] = sort(w_max_idx);
data = unsorted_whisk_heat(:,w_idx_ix)'; 
%set unsampled heatmap to nan; 
[nr,nc] = size(data);
pcolor([data nan(nr,1); nan(1,nc+1)]);
caxis([0 1])
shading flat;
colormap turbo
if strcmp(hilbertVar,'pole')
    set(gca,'xdir','reverse','xtick',1:10:length(touch_x),'xlim',[1 length(touch_x)],'xticklabel',-1:1:1,'ydir','reverse')
elseif strcmp(hilbertVar,'phase')
    set(gca,'xtick',1:10:length(touch_x),'xlim',[1 length(touch_x)],'xticklabel',-1:1:1,'ydir','reverse')
elseif strcmp(hilbertVar,'angle')
    set(gca,'xtick',10:20:length(whisk_x),'xticklabel',-20:20:80,'xlim',[0 length(whisk_x)])
end
title('ix whisk sorted')

%plotting all whisk+touch object location tuned units touch responses sorted by whisk peak 
subplot(3,2,5)
% unsorted_heat_touch = normalize_var(cell2mat(touch_heat(touch_ix_idx)')',0,1);
unsorted_heat_touch = norm_new(cell2mat(touch_heat(touch_ix_idx)')');
data = unsorted_heat_touch(:,w_idx_ix)'; 
%set unsampled heatmap to nan; 
[nr,nc] = size(data);
pcolor([data nan(nr,1); nan(1,nc+1)]);
caxis([0 1])
shading flat;
colormap turbo
if strcmp(hilbertVar,'pole')
    set(gca,'xdir','reverse','xtick',1:10:length(touch_x),'xlim',[1 length(touch_x)],'xticklabel',-1:1:1,'ydir','reverse')
elseif strcmp(hilbertVar,'phase')
    set(gca,'xtick',1:10:length(touch_x),'xlim',[1 length(touch_x)],'xticklabel',-1:1:1,'ydir','reverse')
elseif strcmp(hilbertVar,'angle')
    set(gca,'xtick',10:20:length(whisk_x),'xticklabel',-20:20:80,'xlim',[0 length(whisk_x)])
end
title('ix touch whisk sorted by whisk')

%plotting all whisk+touch object location tuned units whisk responses sorted by touch peak 
subplot(3,2,6)
% unsorted_whisk_heat = normalize_var(cell2mat(whisk_heat(whisk_ix_idx)')',0,1);
unsorted_whisk_heat = norm_new(cell2mat(whisk_heat(whisk_ix_idx)')');
data = unsorted_whisk_heat(:,t_idx_ix)'; 
%set unsampled heatmap to nan; 
[nr,nc] = size(data);
pcolor([data nan(nr,1); nan(1,nc+1)]);
caxis([0 1])
shading flat;
colormap turbo
if strcmp(hilbertVar,'pole')
    set(gca,'xdir','reverse','xtick',1:10:length(touch_x),'xlim',[1 length(touch_x)],'xticklabel',-1:1:1,'ydir','reverse')
elseif strcmp(hilbertVar,'phase')
    set(gca,'xtick',1:10:length(touch_x),'xlim',[1 length(touch_x)],'xticklabel',-1:1:1,'ydir','reverse')
elseif strcmp(hilbertVar,'angle')
    set(gca,'xtick',10:20:length(whisk_x),'xticklabel',-20:20:80,'xlim',[0 length(whisk_x)])
end
title('ix whisk sorted by touch')


    