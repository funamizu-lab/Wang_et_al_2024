function Figure4B_S3B_Rand_trial_activity_20240417


%repeat
process_Rand_trial_activity_20240417('repeat_OFC_20240427');
process_Rand_trial_activity_20240417('repeat_Hippo_20240427');
process_Rand_trial_activity_20240417('repeat_AC_20240427');
process_Rand_trial_activity_20240417('repeat_PPC_20240427');

%zigzag
process_Rand_trial_activity_20240417('altern_OFC_20240427');
process_Rand_trial_activity_20240417('altern_Hippo_20240427');
process_Rand_trial_activity_20240417('altern_AC_20240427');
process_Rand_trial_activity_20240417('altern_PPC_20240427');
process_Rand_trial_activity_20240417('altern_M1_20240427');
process_Rand_trial_activity_20240417('altern_STR_20240427');

end



function process_Rand_trial_activity_20240417(folders)

tic

[analysis_dir,depth_def] = eval(folders);


shuff_Number = 100;
all_ave_data = cell(length(analysis_dir),1);
all_rand_1half = cell(length(analysis_dir), shuff_Number);
all_rand_2half = cell(length(analysis_dir), shuff_Number);

parfor i = 1:length(analysis_dir)
    disp([i,length(analysis_dir)]);
    temp_dir = analysis_dir{i};
    [~,average_data,rand_average_data_1half,rand_average_data_2half] = getData(temp_dir,depth_def,shuff_Number);
    if (~isempty(average_data))
        all_ave_data{i} = average_data;
        for j = 1 : shuff_Number
            all_rand_1half{i,j} = rand_average_data_1half{j};
            all_rand_2half{i,j}  = rand_average_data_2half{j};
        end
    end
end

ept = any(~cellfun('isempty',all_rand_1half),2); 

all_ave_data = all_ave_data(ept,:);
all_rand_1half = all_rand_1half(ept,:);
all_rand_2half = all_rand_2half(ept,:);



all_average_data=cell2mat(all_ave_data);
[~,maxOrder] =max(all_average_data,[],2);
[~,sortOrder]=sort(maxOrder);
all_average_data = all_average_data(sortOrder,:);

all_rand_average_data_1half=cell(shuff_Number,1);
all_rand_average_data_2half=cell(shuff_Number,1);

for i=1:shuff_Number
    temp1 = cell2mat(all_rand_1half(:,i));
    all_rand_average_data_1half{i}=temp1(sortOrder,:);
    temp2 = cell2mat(all_rand_2half(:,i));
    all_rand_average_data_2half{i}=temp2(sortOrder,:);
end
toc

Regress_coeff = RandHeapmat(all_average_data, all_rand_average_data_1half,all_rand_average_data_2half,shuff_Number);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [frameSound,raw_average_data,rand_1half_data,rand_2half_data] = getData(pathname,depth_def,shuff_number)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(pathname)
temp = dir('Supple_20240415_norm_trace*');
if length(temp) ~= 1
    disp(temp)
    hoge
end
load(temp.name);
%norm_spike max_time spike01 mean_spike std_spike

temp = dir('sig_HMM_neurons_20230310*');
if length(temp) ~= 1
    disp(temp)
    hoge
end
load(temp.name);

temp = dir('depth_spike_20230427*');
if length(temp) == 1
    load(temp.name);
   
    if size(p_task,1) ~= length_neuron
        disp([length(p_task) length_neuron])
        hoge
    end
    
    if depth_def == 1
        depth_neuron = find(spike_depth <= def_depth(1));
    else
        depth_neuron = find(spike_depth > def_depth(1) & spike_depth <= def_depth(2));
    end
elseif isempty(temp)
    depth_neuron = 1:size(p_task,1); %Use all the neurons
    hoge
else
    disp(temp)
    hoge
end

frameSound = frame.frame_sound;

threshold = 10;
sig_sound = get_sig_neuron_all(p_task,threshold);
sig_choice = get_sig_neuron_all(p_task2,threshold);
sig_neuron_all = union(sig_sound,sig_choice);
sig_neuron_all = intersect(sig_neuron_all,depth_neuron);

use_sig_neuron = sig_neuron_all;

if (~isempty(use_sig_neuron))
    Data= get_half_spike(neuron,shuff_number,use_sig_neuron);
    grand_average= Data.raw_average_data;
    temp_1half=Data.shuff_1half;
    temp_2half=Data.shuff_2half;
else
    grand_average = [];
    temp_1half = [];
    temp_2half = [];
end

rand_aveData_1half= temp_1half;
rand_aveData_2half = temp_2half;


grand_ave_data=grand_average;
[~,maxOrder] =max(grand_ave_data,[],2);
[~,sortOrder]=sort(maxOrder);
raw_average_data = grand_ave_data(sortOrder,:);

if (~isempty(rand_aveData_1half) && ~isempty(rand_aveData_2half))
    rand_1half_data=cell(shuff_number,1);
    rand_2half_data=cell(shuff_number,1);
    for i=1:shuff_number
        temp1 = cell2mat(rand_aveData_1half(:,i));
        rand_1half_data{i}=temp1(sortOrder,:);
        temp2 = cell2mat(rand_aveData_2half(:,i));
        rand_2half_data{i}=temp2(sortOrder,:);
    end
else
    raw_average_data = [];
    rand_1half_data = [];
    rand_2half_data = [];
end

end


function NormTrace = get_half_spike(neuron,shuff_Number,sig_neuron)


for i=1:length(neuron)
    if(~isempty(neuron(i).norm_spike_filter))
    Trial_number= size(neuron(i).norm_spike_filter,1);
    end
end

shuffTrial = cell(shuff_Number,1);
for j=1:shuff_Number
    shuffTrial{j}=randperm(Trial_number,round(Trial_number/2));
end

grand_average_data = cell(length(sig_neuron),1);
random_1half_average= cell(length(sig_neuron),shuff_Number);
random_2half_average= cell(length(sig_neuron),shuff_Number);
for i=1:length(sig_neuron)
    spike_trace= neuron(sig_neuron(i)).norm_spike_filter;
    tmp = mean(spike_trace,1);
    if max(tmp) ~= 0
        grand_average_data{i} = rescale(tmp);
    else
        grand_average_data{i} = tmp;
    end
    for n=1:shuff_Number
        tmp=setdiff(1:Trial_number,shuffTrial{n});
        random_1half_average{i,n} = rescale(mean(spike_trace(shuffTrial{n},:),1));
        random_2half_average{i,n} = rescale(mean(spike_trace(tmp,:),1));
    end
end
randshuff_1half = cell2mat(random_1half_average);
randshuff_2half = cell2mat(random_2half_average);


time = size(grand_average_data{1},2);

neuron_num = length(grand_average_data);
NormTrace.raw_average_data = cell2mat(grand_average_data);
NormTrace.shuff_1half= mat2cell(randshuff_1half,neuron_num,time*ones(1,shuff_Number));
NormTrace.shuff_2half= mat2cell(randshuff_2half,neuron_num,time*ones(1,shuff_Number));
end

function sig_neuron = get_sig_neuron_all(p_task,new_p_thre)

use_p_sound = min(p_task,[],2);
use_p_sound = -log10(use_p_sound);
sig_neuron = find(use_p_sound > new_p_thre);

end


function Regress_coeff=RandHeapmat(all_average_data, all_rand_average_data_1half,all_rand_average_data_2half,shuff_Number)


frame_sound = 2000;
% soundPeriod = frame_sound+(-99:1100);
allPeriod = frame_sound+(-1999:4000);

usePeriod = allPeriod;
raw_data=all_average_data(:,usePeriod);
for i=1:size(raw_data,1)
if(max(raw_data(i,:))~=min(raw_data(i,:)))
    raw_data(i,:)=rescale(raw_data(i,:));
else    
    raw_data(i,:)=zeros(1,length(usePeriod));
end
end
[~,maxOrder] =max(raw_data,[],2);
[~,sortOrder]=sort(maxOrder);
raw_data = raw_data(sortOrder,:);
% [~,x_sound]=max(raw_data,[],2);
y_sound=1:size(raw_data,1);

Regress_coeff= zeros(shuff_Number,1);
x_sound_test1 = cell(shuff_Number,2);
x_sound_test2 = cell(shuff_Number,2);
for i=1:shuff_Number
    
    test1=all_rand_average_data_1half{i};

    test_data1=test1(:,usePeriod);
    test_data1=test_data1(sortOrder,:);
    [~,max_1half]=max(test_data1(y_sound,:),[],2);
    x_sound_test1{i,2} =max_1half';
    
    test2=all_rand_average_data_2half{i};

    test_data2=test2(:,usePeriod);
    test_data2=test_data2(sortOrder,:);
    [~,max_2half]=max(test_data2(y_sound,:),[],2);
    x_sound_test2{i,2} =max_2half';
    
    Regress_coeff(i)=corr(max_1half,max_2half,'Type','Spearman');
end



figure
imagesc(raw_data);
hold on
xlim([0,size(raw_data,2)]);
ylim([1,size(raw_data,1)]);
yticks([1,size(raw_data,1)]);
xticks([0,2000, size(raw_data,2)]);
xticklabels({'-2000','0','4000'});
xlabel('Time from sound onset');
ylabel('Neuron number');
title('raw data');
set(gca, 'YDir','reverse')
set(gcf,'Position',[250 250 215 208])


data2=cell2mat(x_sound_test1(:,2));

h=figure('Position',[250 250 215 208]);
errorbar(y_sound,mean(data2),std(data2),'CapSize',0,'Color',[.5 .5 .5]);
hold on
scatter(y_sound,mean(data2),5,'k','filled');
hold on
% plot(y_sound,x_sound,'r');
ylim([0 6000]);
yticks([0 2000 6000]);
yticklabels({'-2000','0','4000'});
xlim([1 size(data2,2)]);
xticks([1 size(data2,2)]);
title({['r = ', num2str(round(mean(Regress_coeff),3))]});
xlabel('Neuron number');
ylabel('Time from sound onset');
view(90,90);

set(h,'PaperPositionMode','auto');

end

