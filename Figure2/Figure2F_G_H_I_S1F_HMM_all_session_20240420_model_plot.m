%{
----------------------------------------------------------------------------
Analyzing behavioral data
At least for the correct rate
----------------------------------------------------------------------------
%}

function Figure2F_G_H_I_S1F_HMM_all_session_20240420_model_plot(root_dir)

% root_dir = '/Volumes/Extreme SSD';

repeat_all = HMM_all_session_analysis_20240420_model_plot_repeat(root_dir);
zigzag_all = HMM_all_session_analysis_20240420_model_plot_zigzag(root_dir);


sabun_repeat = repeat_all(:,3) - repeat_all(:,1);
sabun_zigzag = zigzag_all(:,3) - zigzag_all(:,1);

disp('ranksum repeat and alternate')
ranksum(sabun_repeat, sabun_zigzag)

length_repeat = length(sabun_repeat);
length_zigzag = length(sabun_zigzag);
x_repeat = (rand(length_repeat,1) - 0.5) * 0.2 + 1;
x_zigzag = (rand(length_zigzag,1) - 0.5) * 0.2 + 2;

figure
boxplot([sabun_repeat; sabun_zigzag], [ones(length_repeat,1); ones(length_zigzag,1)*2])
hold on
plot(x_repeat, sabun_repeat,'k.')
hold on
plot(x_zigzag, sabun_zigzag,'k.')
hold on
yline(0,':k')
ylabel('△BIC (Hybrid - RL model')


function model_all = HMM_all_session_analysis_20240420_model_plot_repeat(root_dir)

analysis_dir{1} = [root_dir '/Zigzag_repeat_Ephys/Behavior_during_learning/Repeat_first/a16/repeat'];
analysis_dir{2} = [root_dir '/Zigzag_repeat_Ephys/Behavior_during_learning/Repeat_first/W01/repeat'];
analysis_dir{3} = [root_dir '/Zigzag_repeat_Ephys/Behavior_during_learning/Repeat_first/W02/repeat'];
analysis_dir{4} = [root_dir '/Zigzag_repeat_Ephys/Behavior_during_learning/Repeat_first/W08/repeat'];
analysis_dir{5} = [root_dir '/Zigzag_repeat_Ephys/Behavior_during_learning/Repeat_first/W13/repeat'];
analysis_dir{6} = [root_dir '/Zigzag_repeat_Ephys/Behavior_during_learning/Repeat_first/W14/repeat'];
analysis_dir{7} = [root_dir '/Zigzag_repeat_Ephys/Behavior_during_learning/Repeat_first/W22/repeat'];
analysis_dir{8} = [root_dir '/Zigzag_repeat_Ephys/Behavior_during_learning/Repeat_first/W31/repeat'];


[RL, RL_session, ~] = get_BIC_data('HMM_RL_20211029_indivi_session.mat', analysis_dir, 0);
[HMM, HMM_session, ~] = get_BIC_data('HMM_STATE_20211031_indivi_session.mat', analysis_dir, 0);
[hybrid, hybrid_session, ~] = get_BIC_data('HMM_hybrid_20240414_indivi_session.mat', analysis_dir, 1);

memory = [
375.3515992
248.9792507
333.5463673
351.8277441
310.1965652
353.3866807
337.4002101
404.9674936
];

memory_no = [
381.5225825
284.1276891
360.4233116
373.6052593
366.2921693
378.1698069
357.9437678
428.951468
];

model_all = [RL,HMM,hybrid,memory_no,memory];
model_sabun(:,1) = model_all(:,2)-model_all(:,1);
model_sabun(:,2) = model_all(:,3)-model_all(:,1);
model_sabun(:,3) = model_all(:,4)-model_all(:,1);
model_sabun(:,4) = model_all(:,5)-model_all(:,1);

figure
boxplot(model_sabun,'Labels',{'State','Hybrid','Memory','F-Memory'},'Colors', [0 0.5 0])
for i = 1:4
    rand_x = 0.2 * (rand(length(model_sabun(:,1)),1) - 0.5);
    hold on
    plot(rand_x + i, model_sabun(:,i),'k.')
end
ylabel('△BIC from RL model')
yline(0,':k')
ylim([-25,60])
title('Repeat')


RL_HMM = signrank(RL, HMM);
RL_hybrid = signrank(RL, hybrid);
RL_memory_no = signrank(RL, memory_no);
RL_memory = signrank(RL, memory);

p = [RL_HMM, RL_hybrid,RL_memory_no,RL_memory];
disp(p)

% HMM_memory_no = signrank(HMM, memory_no);
% hybrid_memory_no = signrank(hybrid, memory_no);
% memory_memory_no = signrank(memory, memory_no);
% hybrid_memory = signrank(hybrid, memory);

disp('HMM RL')
make_sabun_regression(HMM_session - RL_session, analysis_dir, [-60 60],'g')
xlabel('Number of sessions')
ylabel('△BIC (State-RL)')
hold on
yline(0,':k')
title('Repeat')

disp('hybrid RL')
make_sabun_regression(hybrid_session - RL_session, analysis_dir, [-60 60],'g')
xlabel('Number of sessions')
ylabel('△BIC (Hybrid-RL)')
hold on
yline(0,':k')
title('Repeat')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function make_sabun_regression(sabun1, analysis_dir, c_axis, color)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%sabun1 = HMM_session - RL_session;
b = nan(2,length(analysis_dir));
temp_session = (1:7)';
for i = 1:length(analysis_dir)
    b(:,i) = regress(sabun1(i,:)',[temp_session,ones(7,1)]);
end
use_b = median(b,2);
signrank(b(1,:))

%at 1
use_y(1) = use_b(1) + use_b(2);
use_y(2) = use_b(1)*7 + use_b(2);

figure
plot(sabun1','k')
hold on
plot([1 7], use_y , color)
set(gca,'xlim',[0.5 7.5])
set(gca,'ylim',c_axis)
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [min_BIC, BIC_session, para_set] = get_BIC_data(save_file, analysis_dir, para_analyze)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

min_BIC = nan(length(analysis_dir),1);
BIC_session = nan(length(analysis_dir),7);
para_set = nan(length(analysis_dir),7);

for i = 1:length(analysis_dir)
    cd(analysis_dir{i})
    load(save_file)

    ave_BIC = mean(BIC_all);
    min_BIC(i,1) = min(ave_BIC);
    
    temp = find(ave_BIC == min(ave_BIC));
    %BIC_session(i,:) = min(BIC_all');
    BIC_session(i,:) = BIC_all(:,temp)';
    
    if para_analyze
        for j = 1:7
            temp_para = para_max(j).matrix;
            para_set(i,j) = temp_para(temp,5); %ratio parameters
        end
    end
end

if para_analyze == 0
    para_set = [];
end
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model_all = HMM_all_session_analysis_20240420_model_plot_zigzag(root_dir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

analysis_dir{1}  = [root_dir '/Zigzag_repeat_Ephys/Behavior_during_learning/Zigzag_first/a17/zigzag'];
analysis_dir{2}  = [root_dir '/Zigzag_repeat_Ephys/Behavior_during_learning/Zigzag_first/W03/zigzag'];
analysis_dir{3}  = [root_dir '/Zigzag_repeat_Ephys/Behavior_during_learning/Zigzag_first/W04/zigzag'];
analysis_dir{4}  = [root_dir '/Zigzag_repeat_Ephys/Behavior_during_learning/Zigzag_first/W05/zigzag'];
analysis_dir{5}  = [root_dir '/Zigzag_repeat_Ephys/Behavior_during_learning/Zigzag_first/W06/zigzag'];
analysis_dir{6}  = [root_dir '/Zigzag_repeat_Ephys/Behavior_during_learning/Zigzag_first/W21/zigzag'];
analysis_dir{7}  = [root_dir '/Zigzag_repeat_Ephys/Behavior_during_learning/Zigzag_first/W23/zigzag'];
analysis_dir{8}  = [root_dir '/Zigzag_repeat_Ephys/Behavior_during_learning/Zigzag_first/W26/zigzag'];
analysis_dir{9}  = [root_dir '/Zigzag_repeat_Ephys/Behavior_during_learning/Zigzag_first/W27/zigzag'];
analysis_dir{10} = [root_dir '/Zigzag_repeat_Ephys/Behavior_during_learning/Zigzag_first/W29/zigzag'];
analysis_dir{11} = [root_dir '/Zigzag_repeat_Ephys/Behavior_during_learning/Zigzag_first/W30/zigzag'];

[RL, RL_session, ~] = get_BIC_data('HMM_RL_20211029_indivi_session.mat', analysis_dir, 0);
[HMM, HMM_session, ~] = get_BIC_data('HMM_STATE_20211031_indivi_session.mat', analysis_dir, 0);
[hybrid, hybrid_session, para_hybrid] = get_BIC_data('HMM_hybrid_20240414_indivi_session.mat', analysis_dir, 1);

memory = [
224.3772304
261.1057393
271.4963807
359.0533899
306.5805575
384.160263
377.1150781
382.1546935
353.9521604
415.8044781
365.0163476
];

memory_no = [
242.7625713
283.5016458
275.8954906
375.9561571
318.7735443
392.2144985
384.1909122
397.6752438
381.6709621
425.0103727
370.8888553
];

model_all = [RL,HMM,hybrid,memory_no,memory];
model_sabun(:,1) = model_all(:,2)-model_all(:,1);
model_sabun(:,2) = model_all(:,3)-model_all(:,1);
model_sabun(:,3) = model_all(:,4)-model_all(:,1);
model_sabun(:,4) = model_all(:,5)-model_all(:,1);

figure
boxplot(model_sabun,'Labels',{'State','Hybrid','Memory','F-Memory'},'Colors', 'm')
for i = 1:4
    rand_x = 0.2 * (rand(length(model_sabun(:,1)),1) - 0.5);
    hold on
    plot(rand_x + i, model_sabun(:,i),'k.')
end
ylabel('△BIC from RL model')
hold on
yline(0,':k')
ylim([-25,60])
title('Alternate')

% RL_HMM = signrank(RL, HMM);
% RL_hybrid = signrank(RL, hybrid);
% RL_memory_no = signrank(RL, memory_no);
% RL_memory = signrank(RL, memory);
% p = [RL_HMM, RL_hybrid,RL_memory_no,RL_memory];
% disp(p)

% HMM_memory_no = signrank(HMM, memory_no);
hybrid_memory_no = signrank(hybrid, memory_no);
memory_memory_no = signrank(memory, memory_no);
% hybrid_memory = signrank(hybrid, memory);
disp(hybrid_memory_no)
disp(memory_memory_no)

disp('HMM RL')
make_sabun_regression(HMM_session - RL_session, analysis_dir, [-60 60],'m')
ylabel('△BIC (State - RL model')
xlabel('Number of sessions')
hold on
yline(0,':k')
title('Alternate')

disp('hybrid RL')
make_sabun_regression(hybrid_session - RL_session, analysis_dir, [-60 60],'m')
ylabel('△BIC (Hybrid - RL model')
xlabel('Number of sessions')
hold on
yline(0,':k')
title('Alternate')


make_sabun_regression(para_hybrid, analysis_dir, [0 1],'m')
ylabel('Pmix (hybrid model)')
xlabel('Number of sessions')
title('Alternate')


