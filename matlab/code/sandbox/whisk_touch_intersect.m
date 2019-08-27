% selectedCells = find(cellfun(@(x) isfield(x.meta.touchProperties,'responseWindow'),U)~=0);
selectedCells = cellfun(@(x) strcmp(x.meta.touchProperties.responseType,'excited'),U);
tStruct = object_location_quantification(U,selectedCells,'pole');

wStruct = whisking_location_quantification(U,selectedCells,'pole');
%% intesrection of whisking and touch 

touch_OL = cellfun(@(x) x.is_tuned==1,tStruct);
rc = numSubplots(sum(touch_OL)); 

sel_tstructs = tStruct(touch_OL);
sel_wstructs = wStruct(touch_OL);


figure(99);clf
whisk_touch_pair = cell(1,sum(touch_OL)); 
for g = 1:sum(touch_OL)
    curr_w = sel_wstructs{g}.stim_response.values;
    curr_t = sel_tstructs{g}.stim_response.values;
    
    whisk_x = round(round(min(curr_w(:,1)),1):.1:round(max(curr_w(:,1)),1),1);
    whisk_response = interp1(curr_w(:,1),curr_w(:,2),whisk_x);
    whisk_std = interp1(curr_w(:,1),curr_w(:,3),whisk_x);
    
    touch_x = round(round(min(curr_t(:,1)),1):.1:round(max(curr_t(:,1)),1),1);
    touch_response = interp1(curr_t(:,1),curr_t(:,2),touch_x);
    touch_std = interp1(curr_t(:,1),curr_t(:,3),touch_x);
    
    [~,~,whisk_idx] = intersect(touch_x,whisk_x);
    [overlap_x,~,touch_idx] = intersect(whisk_x,touch_x);
    
    figure(99);subplot(rc(1),rc(2),g)
%     shadedErrorBar(overlap_x,whisk_response(whisk_idx),whisk_std(whisk_idx),'c')
%     hold on; shadedErrorBar(overlap_x,touch_response(touch_idx),touch_std(touch_idx),'b')
%     plot(overlap_x,whisk_response(whisk_idx),'c')
%     hold on; plot(overlap_x,touch_response(touch_idx),'b')
    scatter(whisk_response(whisk_idx),touch_response(touch_idx))
    all_responses = [whisk_response(whisk_idx) touch_response(touch_idx)];
    set(gca,'xlim',[min(all_responses) max(all_responses)],'ylim',[min(all_responses) max(all_responses)])
    hold on; plot([min(all_responses) max(all_responses)],[min(all_responses) max(all_responses)],'-.k')
    
    normed_responses = normalize_var(all_responses,0,1);
    whisk_touch_pair{g} = reshape(normed_responses,numel(normed_responses)./2,2); 
    
    %calculation of difference to quantify the effect touch has on whisking
%     response_difference = touch_response(touch_idx) - whisk_response(whisk_idx);
%     
%     figure(8);subplot(rc(1),rc(2),g)
%     scatter(touch_response(touch_idx),response_difference,'k')

end

all_values = cell2mat(whisk_touch_pair');
figure(12);clf
scatter(all_values(:,1),all_values(:,2),'k')
hold on; plot([0 1],[0 1],'-.k')
set(gca,'xlim',[0 1],'ylim',[0 1])
ylabel('normalized touch responses')
xlabel('normalized whisk responses')
axis square

lm = fitlm(all_values(:,1),all_values(:,2))
    
