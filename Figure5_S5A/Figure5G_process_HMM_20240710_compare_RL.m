
function Figure5G_process_HMM_20240710_compare_RL


process_HMM_20240710_compare_RL('altern_AC_20240427')

return

function process_HMM_20240710_compare_RL(folders, kaiseki_number)

switch nargin
    case 0
        hoge
    case 1
        kaiseki_number = 1;
    case 2
        disp('OK to analyze')
    otherwise
        hoge
end

[analysis_dir,depth_def] = eval(folders);

brainarea = folders(1:10);

%This values are from 'HMM_all_session_analysis_20240423_ephys_integ3'
repeat_min_para_hybrid = [6, 6, 6, 4, 4];
repeat_min_para_hybrid = repeat_min_para_hybrid - 1; %fit to [2:8]
%a17, w03, w06, w22, w31
zigzag_min_para_hybrid = [3, 3, 5, 3, 2, 3, 3, 3, 2];
zigzag_min_para_hybrid = zigzag_min_para_hybrid - 1; %fit to [2:8]
%w02, w13, w14, w23, w25, w27, w28, w29, w30
all_min_hybrid = [repeat_min_para_hybrid, zigzag_min_para_hybrid];

%Based on the directory name, need to pick up the right values
clear BIC_all2 para_max2
length(analysis_dir)
for i = 1:length(analysis_dir)
    %depending on the name of directory, set the parameter number
    temp_dir = analysis_dir{i};
    clear temp
    temp(1) = contains(temp_dir,'a17');
    temp(2) = contains(temp_dir,'W03');
    temp(3) = contains(temp_dir,'W06');
    temp(4) = contains(temp_dir,'W22');
    temp(5) = contains(temp_dir,'W31');
    
    temp(6) = contains(temp_dir,'W02');
    temp(7) = contains(temp_dir,'W13');
    temp(8) = contains(temp_dir,'W14');
    temp(9) = contains(temp_dir,'W23');
    temp(10) = contains(temp_dir,'W25');
    temp(11) = contains(temp_dir,'W27');
    temp(12) = contains(temp_dir,'W28');
    temp(13) = contains(temp_dir,'W29');
    temp(14) = contains(temp_dir,'W30');
    
    if sum(temp) ~= 1
        disp(temp)
        hoge
    else
        temp = temp == 1;
        hybrid_para(i) = all_min_hybrid(temp);
    end
end

%Separate the sessions with the behavioral performance
[~, ~, ~, temp_para1] = get_BIC_data_ephys_20240523('State20240414*.mat', 2:8, analysis_dir, hybrid_para);
para_hybrid = temp_para1(:,5); %size = [31,5]

thre = median(para_hybrid);
disp(thre)

%hybrid_L = (1-para(5)) * relative_Q(1) + para(5) * prob_LH(i,1);
%para is high, use more state-based model 
temp1 = find(para_hybrid >= thre);
temp2 = find(para_hybrid < thre);

%
[all_pre_sound1, right_sig_neuron1,left_sig_neuron1,~] = ...
    session_separete_data(temp1, analysis_dir, kaiseki_number,depth_def,brainarea);

[all_pre_sound2, right_sig_neuron2,left_sig_neuron2,~] = ...
    session_separete_data(temp2, analysis_dir, kaiseki_number,depth_def,brainarea);

disp([length(temp1), length(temp2)])
%para_hybrid


%%%%%%%%%%%%
%Compute the BIC
sig_neuron1 = union(right_sig_neuron1, left_sig_neuron1);
sig_neuron2 = union(right_sig_neuron2, left_sig_neuron2);
correct1 = all_pre_sound1(sig_neuron1,6);
error1 = -all_pre_sound1(sig_neuron1,7);

correct2 = all_pre_sound2(sig_neuron2,6);
error2 = -all_pre_sound2(sig_neuron2,7);

all_correct = [correct1; correct2];
all_error = [error1; error2];
n_data = length(all_correct);

%regression
%[b,~,~,~,stats] = regress(index2,[index1,ones(length(index1),1)]);
[b_regress(1),~,stats] = glmfit(correct1,error1,'normal','link','identity','Constant','off');
p_regress(1) = stats.p;
[b_regress(2),~,stats] = glmfit(correct2,error2,'normal','link','identity','Constant','off');
p_regress(2) = stats.p;

[b_regress(3),~,stats] = glmfit(all_correct,all_error,'normal','link','identity','Constant','off');
p_regress(3) = stats.p;

figure
plot([-1 1],[0 0],'k:')
hold on
plot([0 0],[-1 1],'k:')
hold on
plot([-1 1],[-1 1],'k:')
hold on
plot(correct1, error1, 'm.')
hold on
plot([-1,1],[-b_regress(1),b_regress(1)],'m')
hold on
plot([-1,1],[-b_regress(3),b_regress(3)],Color=[.7 .7 .7])
hold on
plot(correct2, error2, 'k.')
hold on
plot([-1,1],[-b_regress(2),b_regress(2)],'c')
hold on
plot([-1,1],[-b_regress(3),b_regress(3)],'k')

set(gca,'xlim',[-1, 1],'ylim',[-1, 1],'fontname','Arial')
set(gcf,'Position',[584,652,295,263])

title('Alternating condition (AC)')

%Compute the BICs based on the root mean squared error
r_RSS = (b_regress(1) * correct1 - error1) .* (b_regress(1) * correct1 - error1);
%r_RSS = sqrt(sum(r_RSS));

z_RSS = (b_regress(2) * correct2 - error2) .* (b_regress(2) * correct2 - error2);
%z_RSS = sqrt(sum(z_RSS));

rz_RSS = [r_RSS; z_RSS];
rz_RSS = sqrt(sum(rz_RSS));


all_RSS = (b_regress(3) * all_correct - all_error) .* (b_regress(3) * all_correct - all_error);
all_RSS = sqrt(sum(all_RSS));

BIC_rz = n_data * log(rz_RSS/n_data) + 2 * log(n_data);
BIC_all = n_data * log(all_RSS/n_data) + log(n_data);
disp([BIC_rz, BIC_all])

disp([length(correct1), length(correct2)])

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [all_pre_sound, right_sig_neuron,left_sig_neuron,non_sig_neuron] = ...
    session_separete_data(thre_session, analysis_dir, kaiseki_number,depth_def,brainarea)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
all_pre_sound = [];
all_p_pre_sound = [];
all_norm_spike = [];
all_p_surprise = [];
all_sound_correct_error = [];
all_p_sound_correct_error = [];
all_neuron_choice_category = [];

for i = 1:length(thre_session)
    disp([i,length(thre_session)])
    temp_dir = analysis_dir{thre_session(i)};
    
    [sound_neuron, p_sound_neuron, norm_spike, p_surprise, ...
    sound_correct_error, p_sound_correct_error, neuron_choice_category] = ...
        HMM_ephys_20240420_surprise_before_sound_depth(temp_dir, kaiseki_number,depth_def);

    all_pre_sound = [all_pre_sound; sound_neuron];
    all_p_pre_sound = [all_p_pre_sound; p_sound_neuron];
    all_norm_spike = [all_norm_spike; norm_spike];
    all_p_surprise = [all_p_surprise; p_surprise];

    all_sound_correct_error = [all_sound_correct_error; sound_correct_error];
    all_p_sound_correct_error = [all_p_sound_correct_error; p_sound_correct_error];

    all_neuron_choice_category = [all_neuron_choice_category; neuron_choice_category];
end
delete(gcp('nocreate'))

disp([length(all_pre_sound),length(all_norm_spike),length(all_p_surprise)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot about low, high neurons

if contains(brainarea, 'repeat')
    right_neuron = find(all_pre_sound(:,6) >= 0);
    left_neuron = find(all_pre_sound(:,6) < 0);
    sig_sound = find(all_p_pre_sound(:,3) < 0.01); %Sig diff. between "PREVIOUS" left and right
else
    left_neuron = find(all_pre_sound(:,6) >= 0);
    right_neuron = find(all_pre_sound(:,6) < 0);
    sig_sound = find(all_p_pre_sound(:,3) < 0.01); %Sig diff. between "PREVIOUS" left and right
end

right_sig_neuron = intersect(right_neuron, sig_sound);
left_sig_neuron = intersect(left_neuron, sig_sound);
non_sig_neuron = setdiff(1:length(all_pre_sound), sig_sound);
check_neuron = length(right_sig_neuron) + length(left_sig_neuron) + length(non_sig_neuron);
if check_neuron ~= length(all_pre_sound)
    hoge
else
    disp([length(left_sig_neuron),length(right_sig_neuron),length(non_sig_neuron)])
end
size(all_sound_correct_error)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [min_BIC, BIC_session,ave_para,para_session] = get_BIC_data_ephys_20240523(save_file, BIC_number, use_session, hybrid_para)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear BIC_all2 para_max2
for j = 1:length(use_session)
    cd(use_session{j});
    temp_file = dir(save_file);

    if length(temp_file) ~= 1
        hoge
    else
        load(temp_file.name)
        %Use the BIC with hybrid_para
        BIC_all = BIC_all(BIC_number);
        para_max = para_max(BIC_number,:);
        
        BIC_all2(j,1) = BIC_all(hybrid_para(j));
        para_max2(j,:) = para_max(hybrid_para(j),:);
    end
end

BIC_session = BIC_all2;
min_BIC = mean(BIC_all2);
para_session = para_max2;
ave_para = mean(para_session);

return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p, median_value,rho,pval,b_regress,p_regress] = plot_index_abs_index_sig(index1, index2, right_neuron,left_neuron,non_sig_neuron)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sig_neuron = union(left_neuron, right_neuron);

figure
subplot(2,2,1)
plot([-1 1],[0 0],'k:')
hold on
plot([0 0],[-1 1],'k:')
hold on
plot([-1 1],[-1 1],'k:')
hold on
%plot(index1, index2, 'b.')
plot(index1(non_sig_neuron), index2(non_sig_neuron), 'k.')
hold on
plot(index1(left_neuron), index2(left_neuron), 'b.')
hold on
plot(index1(right_neuron), index2(right_neuron), 'r.')

[rho(1),pval(1)] = corr(index1,index2);
[rho(2),pval(2)] = corr(index1(sig_neuron),index2(sig_neuron));

%regression

[b_regress(1),~,stats] = glmfit(index1,index2,'normal','link','identity','Constant','off');
p_regress(1) = stats.p;
[b_regress(2),~,stats] = glmfit(index1(sig_neuron),index2(sig_neuron),'normal','link','identity','Constant','off');
p_regress(2) = stats.p;

subplot(2,2,2)
plot([0 1],[0 1],'k:')
hold on
plot(abs(index1(non_sig_neuron)), abs(index2(non_sig_neuron)), 'k.')
hold on
plot(abs(index1(left_neuron)), abs(index2(left_neuron)), 'b.')
hold on
plot(abs(index1(right_neuron)), abs(index2(right_neuron)), 'r.')

p = signrank(abs(index1), abs(index2));
median_value = [median(abs(index1)), median(abs(index2))];

subplot(2,2,3)
plot([-1 1],[0 0],'k:')
hold on
plot([0 0],[-1 1],'k:')
hold on
plot([-1 1],[-1 1],'k:')
hold on
plot(index1(left_neuron), index2(left_neuron), 'b.')
hold on
plot(index1(right_neuron), index2(right_neuron), 'r.')
subplot(2,2,4)
plot([0 1],[0 1],'k:')
hold on
plot(abs(index1(left_neuron)), abs(index2(left_neuron)), 'b.')
hold on
plot(abs(index1(right_neuron)), abs(index2(right_neuron)), 'r.')

return
