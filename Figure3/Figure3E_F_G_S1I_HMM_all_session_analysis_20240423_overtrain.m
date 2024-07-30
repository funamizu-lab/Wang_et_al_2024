%{
----------------------------------------------------------------------------
Analyzing behavioral data
At least for the correct rate
----------------------------------------------------------------------------
%}

function Figure3E_F_G_S1I_HMM_all_session_analysis_20240423_overtrain

disp('start repeat')
[repeat_session, para_repeat, ave_para_repeat, repeat_sabun] = HMM_all_session_analysis_20240423_each('repeat_overtrain_20240420');
disp('start alternate')
[zigzag_session, para_zigzag, ave_para_zigzag, zigzag_sabun] = HMM_all_session_analysis_20240423_each('zigzag_overtrain_20240420');

disp([median(para_repeat), median(para_zigzag)])
disp([median(ave_para_repeat), median(ave_para_zigzag)])

temp_box = [ones(length(repeat_session),1);2*ones(length(zigzag_session),1)];

figure

boxplot([repeat_sabun(:,1);zigzag_sabun(:,1)],temp_box)
hold on
plot(0.2*(rand(length(repeat_session),1)-0.5) + 1,repeat_sabun(:,1),'k.')
hold on
plot(0.2*(rand(length(zigzag_session),1)-0.5) + 2,zigzag_sabun(:,1),'k.')
hold on
yline(0,':k')
ylim([-40,40])
ylabel('△BIC (HMM - RL model)')
title('Overtrained phase')

figure
boxplot([repeat_sabun(:,2);zigzag_sabun(:,2)],temp_box)
hold on
plot(0.2*(rand(length(repeat_session),1)-0.5) + 1,repeat_sabun(:,2),'k.')
hold on
plot(0.2*(rand(length(zigzag_session),1)-0.5) + 2,zigzag_sabun(:,2),'k.')
ylabel('△BIC (Hybrid - RL model)')
yline(0,':k')
ylim([-50,20])
title('Overtrained phase')

ranksum(repeat_sabun(:,1),zigzag_sabun(:,1))
ranksum(repeat_sabun(:,2),zigzag_sabun(:,2))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mouse_session, para_hybrid, ave_para_hybrid, model_sabun] = HMM_all_session_analysis_20240423_each(task_condition)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

all_directory = eval(task_condition);

RL = nan(length(all_directory),1);
HMM = nan(length(all_directory),1);
hybrid = nan(length(all_directory),1);
memory = nan(length(all_directory),1);
memory_no = nan(length(all_directory),1);
min_para_hybrid = nan(length(all_directory),1);
mouse_session = nan(1,5);
ave_para_hybrid = nan(length(all_directory),5);

RL_session = [];
HMM_session = [];
hybrid_session = [];
memory_session = [];
memory_no_session = [];

para_hybrid = [];
para_memory = [];

for i = 1:length(all_directory)
    disp(i)
    use_session = all_directory(i).matrix;
    mouse_session(i) = length(use_session);
    
    [RL(i,:), temp_RL] = get_BIC_data_ephys('RL20230310*.mat', 1:7, use_session);
    [HMM(i,:), temp_HMM] = get_BIC_data_ephys('State20230310*.mat', 1:7, use_session);
    [hybrid(i,:), temp_hybrid, ave_para_hybrid(i,:), temp_para1,temp_min] = get_BIC_data_ephys('State20240414*.mat', 2:8, use_session);

    [memory(i,:), temp_memory, ~, temp_para2] = get_BIC_data_ephys('GAO20240416*.mat', 1:5, use_session);
    [memory_no(i,:), temp_memory_no] = get_BIC_data_ephys('GAO20240418_nofor*.mat', 1:5, use_session);

    RL_session = [RL_session; temp_RL];
    HMM_session = [HMM_session; temp_HMM];
    hybrid_session = [hybrid_session; temp_hybrid];
    memory_session = [memory_session; temp_memory];
    memory_no_session = [memory_no_session; temp_memory_no];
    
    para_hybrid = [para_hybrid; temp_para1];
    para_memory = [para_memory; temp_para2];
    
    min_para_hybrid(i,1) = temp_min;
end

session_all = [];
for i = 1:length(mouse_session)
    temp = mouse_session(i);
    session_all = [session_all; ones(temp,1)*i];
end

model_all = [RL,HMM,hybrid,memory_no,memory];
model_sabun(:,1) = model_all(:,2)-model_all(:,1);
model_sabun(:,2) = model_all(:,3)-model_all(:,1);
model_sabun(:,3) = model_all(:,4)-model_all(:,1);
model_sabun(:,4) = model_all(:,5)-model_all(:,1);

% mean(model_all)

ave_para_hybrid = ave_para_hybrid(:,5);
para_hybrid = para_hybrid(:,5);

if contains(task_condition, 'repeat')
    tempcolor = [0 0.5 0];
    tmptitle = 'Repeat';
else
    tempcolor = 'm';
    tmptitle = 'Alterante';
end

figure
boxplot(model_sabun, 'Labels',{'State','Hybrid','Memory','F-Memory'},'Colors', tempcolor)
for i = 1:4
    rand_x = 0.2 * (rand(length(model_sabun(:,1)),1) - 0.5);
    hold on
    plot(rand_x + i, model_sabun(:,i),'k.')
end
yline(0,':k')
ylim([-60,75])
ylabel('△BIC from RL model')
title(tmptitle)


p_lme(1,1) = compute_lme_p(RL_session-HMM_session, session_all);
p_lme(2,1) = compute_lme_p(RL_session-hybrid_session, session_all);
p_lme(3,1) = compute_lme_p(RL_session-memory_no_session, session_all);
p_lme(4,1) = compute_lme_p(RL_session-memory_session, session_all);

p_lme(5,1) = compute_lme_p(memory_no_session-HMM_session, session_all);
p_lme(6,1) = compute_lme_p(hybrid_session-memory_no_session, session_all);
p_lme(7,1) = compute_lme_p(memory_session-memory_no_session, session_all);

p_lme(8,1) = compute_lme_p(hybrid_session-HMM_session, session_all);
p_lme(9,1) = compute_lme_p(hybrid_session-memory_session, session_all);
p_lme(10,1) = compute_lme_p(memory_session-HMM_session, session_all);

tmpp = [p_lme(8,:), p_lme(6,:),p_lme(9,:)];

if contains(task_condition, 'repeat')
    disp(p_lme(1:4))
else
    disp(tmpp)
end

temp = [RL_session;hybrid_session];
min_temp = min(temp);
max_temp = max(temp);


figure
plot([min_temp,max_temp],[min_temp,max_temp],'k:')
hold on
plot(RL_session,hybrid_session,'Color',tempcolor, 'LineStyle','none','Marker','.')
xlabel('BIC of RL model')
ylabel('BIC of hybrid model')
title(tmptitle)

% p = signrank(RL_session, hybrid_session);
% disp(p)


return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = compute_lme_p(sabun_session, session_all)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_name{1} = 'value';
data_name{2} = 'random1';
%linear mixed effect model
data = [sabun_session,session_all];
tbl = table(data(:,1),data(:,2),'VariableNames',data_name);
lme = fitlme(tbl,'value ~ 1+(1|random1)');
p = lme.Coefficients.pValue;
  
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [min_BIC, BIC_session,ave_para,para_session,temp_min] = get_BIC_data_ephys(save_file, BIC_number, use_session)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear BIC_all2 para_max2

for j = 1:length(use_session)
    cd(use_session{j});
    temp_file = dir(save_file);

    if length(temp_file) ~= 1
        hoge
    else
        load(temp_file.name)
        BIC_all2(j,:) = BIC_all(BIC_number);
        para_max2(j).matrix = para_max(BIC_number,:);
    end
end

ave_BIC = mean(BIC_all2);
temp_min = find(ave_BIC == min(ave_BIC),1);

BIC_session = BIC_all2(:,temp_min);
min_BIC = ave_BIC(temp_min);

clear para_session

for j = 1:length(use_session)
    temp_para = para_max2(j).matrix;
    para_session(j,:) = temp_para(temp_min,:);
end
ave_para = mean(para_session);

temp_min = BIC_number(temp_min);
return


