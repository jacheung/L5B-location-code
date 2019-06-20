
touchWindow = [-25:50];
preDecisionTouches = preDecisionTouchMat(U);
phaseBins = 12; %would recommend even numbers;
for rec = 1:length(U)
 [tVar] = atTouch_sorter(U{rec},touchWindow,preDecisionTouches{rec});
 
 phase = tVar.allTouches.S_ctk(:,5);
 if isfield(U{rec}.meta,'responseWindow')
%      tResponse = U{rec}.meta.responseWindow+find(touchWindow==0) ;
     tResponse = [-25 -5] + find(touchWindow==0); %control to find off touch responses
 else
     tResponse = [5 35] + find(touchWindow==0);
 end
  response = mean(tVar.allTouches.R_ntk(:,tResponse(1):tResponse(2)),2);
 
  
  hilbertTouch.S_ctk{rec} = tVar.allTouches.S_ctk(:,[3, 4, 5, 1]);
  hilbertTouch.R_ntk{rec} = response; 
  
end

%%
for k = 3
    phase = hilbertTouch.S_ctk{k}(:,3); 
    amp = hilbertTouch.S_ctk{k}(:,1); 
    midpoint = hilbertTouch.S_ctk{k}(:,2); 
    angle = hilbertTouch.S_ctk{k}(:,4); 
    
    
    responses = normalize_var(hilbertTouch.R_ntk{k},.05,.95);
    [responseSort,idx] = sort(responses); 
    
    featX = phase(idx);
    featY = amp(idx);
    featZ = angle(idx); 
    
    colors = repmat(1 - responseSort,1,3); 
    figure(2332);clf
    
    scatter(featX,featY,100,colors,'filled');
    
end
    