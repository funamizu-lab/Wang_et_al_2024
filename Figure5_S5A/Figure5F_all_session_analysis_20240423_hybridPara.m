%{
----------------------------------------------------------------------------
Analyzing behavioral data
At least for the correct rate
----------------------------------------------------------------------------
%}

function Figure5F_all_session_analysis_20240423_hybridPara

HMM_all_session_analysis_20240423_each('zigzag_overtrain_20240420');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function HMM_all_session_analysis_20240423_each(task_condition)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

all_directory = eval(task_condition);
para_hybrid = [];

for i = 1:length(all_directory)
    disp(i)
    use_session = all_directory(i).matrix;
    [~, ~, ~, temp_para1,~] = get_BIC_data_ephys('State20240414*.mat', 2:8, use_session);    
    para_hybrid = [para_hybrid; temp_para1];   
end

figure
boxplot(para_hybrid(:,5),'Symbol', '', 'Colors', 'm')
hold on
rand_x = 1 + 0.2 * (rand(length(para_hybrid(:,5)),1) - 0.5);
plot(rand_x,para_hybrid(:,5),'k.')
ylim([-0.05,1])
ylabel('Pmix (hybrid model')
title('Alterante')



return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [min_BIC, BIC_session,ave_para,para_session,temp_min] = get_BIC_data_ephys(save_file, BIC_number, use_session)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear BIC_all2 para_max2
BIC_all2 = nan(length(use_session),7);

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
para_session = nan(length(use_session),5);
for j = 1:length(use_session)
    temp_para = para_max2(j).matrix;
    para_session(j,:) = temp_para(temp_min,:);
end
ave_para = mean(para_session);

temp_min = BIC_number(temp_min);
return


