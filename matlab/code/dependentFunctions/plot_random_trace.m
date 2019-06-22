% view traces with whisking, spikes, and touch responses

for k = 1
    currentCell = U{k};
    trialSample = datasample(1:U{k}.k,1);
    whisk = currentCell.S_ctk(1,:,trialSample);
    touchesOn = [find(currentCell.S_ctk(9,:,trialSample)==1) find(currentCell.S_ctk(12,:,trialSample)==1)];
    touchesOff = [find(currentCell.S_ctk(10,:,trialSample)==1) find(currentCell.S_ctk(13,:,trialSample)==1)];
    
    spikes = find(currentCell.R_ntk(1,:,trialSample)==1);
    figure(2310);clf
    if ~isempty(touchesOn) && ~isempty(spikes)
        plot(1:currentCell.t,whisk,'g')
        for g = 1:length(touchesOn)
            tp = touchesOn(g):touchesOff(g); 
        hold on; plot(tp,whisk(tp),'r');
        end
        hold on; scatter(spikes,ones(length(spikes),1)*mean(whisk),'ko','filled')
    end
    
    title(['cell num ' num2str(k) ' at trial num ' num2str(trialSample)])
end
    