% assign a bin number to each trial
clear
tic
filename='../data/validation2.h5';
file_info = h5info(filename);
group_name = file_info.Groups.Name;
tr = h5read(filename, [group_name,'/trials']);
envt = h5read(filename, [group_name,'/trial_info'], [1 1], [2 Inf]);
num_trials=size(tr,2);
difficulties=zeros(1,num_trials);
bin_numbers=difficulties;
acc=zeros(1,h.NumBins);
load('../data/hist_difficulties_validation2_12bins.mat')

for trn=1:num_trials
    [lst,rst,last_cp]=tr{:,trn}; end_state=envt(2,trn);
    if isempty(last_cp)
        time=0;
    else
        time=last_cp(end);
    end 
    num_left=sum(lst>=time); num_right=sum(rst>=time);
    if end_state==1
        diff=num_left-num_right;
    elseif end_state==-1
        diff=num_right-num_left;
    else
        warning('trial with 0 end state')
    end
    difficulties(trn)=diff;
    bin_numbers(trn)=find(diff>=h.BinEdges,1,'last');
    
    % decide with several models...
end


% compute acc on each bin


% compute PP on each bin