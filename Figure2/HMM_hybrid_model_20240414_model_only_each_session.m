%{
----------------------------------------------------------------------------
%The previous 20211029 was not good, maybe the threshold bias moved to the
opposite side
----------------------------------------------------------------------------
%}

function [ave_likeli,BIC_all,log_likeli,para_max,N_trial] = ...
    HMM_hybrid_model_20240414_model_only_each_session(Correct_trial, Chosen_trial, binary_tone, Reward_LCR)

% switch nargin
%     case 0
%         [moto_filename1, pathname1]=uigetfile('*.mat','Block_mat');
%         filename1 = [pathname1, moto_filename1];
%         load(filename1)
%         save_file = ['HMM_State_20211031_',moto_filename1];
%     case 1
%         load(filename1)
%         save_file = ['HMM_State_20211031_',filename1];
%     otherwise
%         hoge
% end
% 
% [Choice_trial,tone_evidence,trial_evidence,use_trial_remove_first,...
%     low,high,correct,error,flip_tone,number_use_trial,number_use_trial_remove_first,...
%     binary_tone,right_trial_all,number_trial_all,right_trial,number_trial] ...
%  = HMM_get_basic_task_structure_20210514(filename1);

%binary_tone is already flipped!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Use all trials for analysis
N_trial = length(Correct_trial);
% Correct_trial = Correct_side(Choice_trial);
% Chosen_trial = Chosen_side(Choice_trial);
% binary_tone = binary_tone(Choice_trial);
% Reward_LCR = Reward_LCR(Choice_trial,[1,3]); %choose chosen trials

clear init_Q
init_Q = Reward_LCR(1,:); %2.4 2.4
    
%reward amount for each trial

reward_trial = zeros(N_trial,1);
for i = 1:N_trial
    %Reward trials
    if Correct_trial(i) == Chosen_trial(i)
        reward_trial(i) = Reward_LCR(i,Chosen_trial(i)+1);
    end
end

%Parameter
bin_x = [0:0.001:1];
stim_LR   = [0,0.25,0.45,0.55,0.75,1]; %3stim
stim_prob_LR   = [0.5,0.25,0.25,0.25,0.25,0.5]; %3stim
stim_prob_kitei = [0.25,0.125,0.125,0.125,0.125,0.25]; %3stim
%Sense_std = [0.01:0.05:1];
%reward = [3,1]; %left-right
%stim_prob = [0.5, 0.5];

%Instead of updating values for each trial, updaiting the state transition
% %Parameter
% para(1) = 0.1; % learning of state transition
% para(2) = 0.2;  % std
% para(3) = 0.3683;  % beta
% para(4) = 0;  % bias
% %para(5) = 0.1;  %inverse temperature for prior

use_para(1).matrix = [2 5];
use_para(2).matrix = [1 2 5];
use_para(3).matrix = [1 2 3 5];
use_para(4).matrix = [1 2 4 5];
use_para(5).matrix = [1 2 3 4 5];

para_temp = [0    0.2 0.5 0 0.1;
             0    0.5 1   0 0.1;
             0    1   0.5 0 0.1;
             0.01 0.2 1   0 0.1;
             0.01 0.5 0.5 0 0.1;
             0.1  0.2 0.5 0 0.1;
             0.1  0.5 1   0 0.1;
             ];
[size_y,size_x] = size(para_temp);
         
%Get maximum prediction
opt = optimset('Display','off');
clear ave_likeli BIC_all log_likeli
for j = 1:length(use_para)
    j
    temp_use_para = use_para(j).matrix;
    parfor i = 1:size_y
    %for i = 1:size_y
        i
        [X1(i,:),FCAL1(i,:),EXITFLAG,OUTPUT] = ...
            fminsearch(@HMM_hybrid_update_211031_all,para_temp(i,:), opt, Correct_trial, reward_trial, Chosen_trial, binary_tone, bin_x, stim_LR, stim_prob_LR, init_Q, temp_use_para);
    end
    min_FCAL = find(FCAL1 == min(FCAL1),1);
    temp_para_max = X1(min_FCAL,:);

    [ave_likeli(1,j),likelihood,Q,Choice_prob,para_max(j,:),state_prob] = ...
        HMM_Hybrid_update_211031_max_all(temp_para_max, Correct_trial, reward_trial, Chosen_trial, binary_tone, bin_x, stim_LR, stim_prob_LR, init_Q, temp_use_para);
    
    sumlog_likeli = sum(log(likelihood));
    log_likeli(1,j) = sumlog_likeli;
    BIC_all(1,j) = -2 * sum(log(likelihood)) + length(temp_use_para) * log(N_trial);
        %log_likelihood(filecount).matrix = log(likelihood);
end
ave_likeli
BIC_all
log_likeli
para_max
N_trial

% figure
% plot(state_prob,'b')
% % hold on
% % plot(Q(:,2),'r')
% set(gca,'xlim',[1 N_trial])
% set(gca,'xtick',[0:40:N_trial])
% set(gca,'ylim',[0 4])
% 
% save_file = 'HMM_all_session_STATE_20211031.mat';
% save(save_file, 'ave_likeli','BIC_all','log_likeli','para_max','N_trial','state_prob')
% 
%hoge

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ave_likeli,likelihood,prob_LH,Choice_prob,para,state_prob] = ...
    HMM_Hybrid_update_211031_max_all(para, Correct_trial, reward, action, Stim_freq, bin_x, stim_LR, stim_prob_LR, init_Q, use_para)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% para(1) = 0.1; % learning rate of state transition
% para(2) = 0.2;  % std
% para(3) = 0.3683;  % beta
% para(4) = 0;  % bias
% %para(5) = 0.1;  %inverse temperature for prior
% use_para(1).matrix = [2];
% use_para(2).matrix = [1 2];
% use_para(3).matrix = [1 2 3];
% use_para(4).matrix = [1 2 4];
% use_para(5).matrix = [1 2 3 4];

temp = find(use_para == 1, 1);
if isempty(temp)
    para(1) = 0;
end
temp = find(use_para == 3, 1);
if isempty(temp)
    softmax_check = 0;
else
    softmax_check = 1;
end
temp = find(use_para == 4, 1);
if isempty(temp)
    para(4) = 0;
end
temp = find(use_para == 5, 1);
if isempty(temp)
    para(5) = 0;
end

%para shusei
if para(1) < 0.001
    para(1) = 0.001;
elseif para(1) > 1
    para(1) = 1;
end
if para(2) < 0
    para(2) = 0;
end
if para(3) < 0
    para(3) = 0;
end
if para(5) < 0
    para(5) = 0;
elseif para(5) > 1
    para(5) = 1;
end
%Plot the gaussian based class for left and right
%Write the figures for bayes computation
N_trial = length(reward);
%long stimulus
for i = 1:length(stim_LR)
    temp_x = normpdf(bin_x,stim_LR(i),para(2));
    temp_x = temp_x ./ sum(temp_x);
    LR_stim(i,:) = stim_prob_LR(i) .* temp_x;
    norm_pdf_long(i,:) = temp_x;
end
%Combination of left and right
L_long = sum(LR_stim([1:3],:));
R_long = sum(LR_stim([4:6],:));
% figure
% plot(L_long,'b')
% hold on
% plot(R_long,'r')

clear Prior Q_left Q_right Choice_prob
%Belief_basis for left choice reward probability
Q = zeros(N_trial,2);
Q(1,:)  = init_Q;
prob_LH = nan(N_trial,2);
prob_LH(1,:)  = [0.5,0.5];
state_prob = nan(N_trial,1);
state_prob(1) = 0.5;
likelihood = zeros(N_trial,1);
Choice_prob = zeros(N_trial,2);
action = action + 1; %0,1 -> 1,2

for i = 1:N_trial
    clear Posterior_LR Decision_LR
    
    %Based on the prob_LH, estimate the posterior probabibility of sound
    %stimulus
    relative_Q = Q(i,:) ./ sum(Q(i,:));
    %Compute hybrid value
    hybrid_L = (1-para(5)) * relative_Q(1) + para(5) * prob_LH(i,1);
    hybrid_R = (1-para(5)) * relative_Q(2) + para(5) * prob_LH(i,2);
    
% %     temp_L = prob_LH(i,1) .* L_long;
% %     temp_R = prob_LH(i,2) .* R_long;
%     temp_L = hybrid_L .* L_long;
%     temp_R = hybrid_R .* R_long;
%     temp_thre = find(temp_R >= temp_L,1);
%     temp_thre = bin_x(temp_thre);
%     if isempty(temp_thre)
%         temp_thre = 1;
%     end
    
    if softmax_check
%         %choice_prior = softmax(prob_LH(i,:),para(3));
%         choice_prior = softmax([temp_thre,1-temp_thre],para(3));
        choice_prior = softmax([hybrid_L,hybrid_R],para(3));
    else
        %choice_prior = temp_thre;
        choice_prior = [hybrid_L,hybrid_R]./sum(hybrid_L+hybrid_R);
    end
    decision_x = choice_prior(1);

    left_choice  = normcdf(decision_x + para(4),Stim_freq(i),para(2)) - normcdf(0,Stim_freq(i),para(2));
    right_choice = normcdf(1,Stim_freq(i),para(2)) - normcdf(decision_x + para(4),Stim_freq(i),para(2));
    left_choice = left_choice / (left_choice + right_choice);
    right_choice = 1 - left_choice;

    if left_choice > 1
        left_choice = 1;
        right_choice = 0;
    elseif right_choice > 1
        left_choice = 0;
        right_choice = 1;
    end
    
    Choice_prob(i,:) = [left_choice,right_choice];
    %Get probability from action and reward;
    likelihood(i) = Choice_prob(i,action(i));
    
    %Forgetting Q learning
    Q(i+1,:) = Qlearning_ori(Q(i,:), action(i), reward(i), para(1));
    
    if i == 1
        prob_LH(i+1,:) = prob_LH(i,:); %No history of updating
        state_prob(i+1,:) = state_prob(i,:); %No history of updating
    else
        pre_state = Correct_trial(i-1);
        post_state = Correct_trial(i);
        
        if pre_state == post_state
            temp_tran = 0;
        else
            temp_tran = 1;
        end
        %Analyzing state transition
        state_prob(i+1,:) = Transition_learning_ori(state_prob(i,:), temp_tran, para(1));
        %Based on the state transition estimate the next stimulus
        if post_state == 1 %High stimulus
            prob_LH(i+1,:) = [state_prob(i+1,:), 1-state_prob(i+1,:)];
        else
            prob_LH(i+1,:) = [1-state_prob(i+1,:), state_prob(i+1,:)];
        end
    end
end

% subplot(2,1,1)
% plot(state_prob,'k')
% subplot(2,1,2)
% plot(prob_LH)
% refreshdata
% drawnow
    
%likelihood keisan
log_likeli = sum(log(likelihood));
log_likeli = log_likeli / N_trial;
ave_likeli = exp(log_likeli);
ave_likeli = -ave_likeli;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ave_likeli_mean = HMM_hybrid_update_211031_all(para, Correct_trial, reward_trial, Chosen_trial, binary_tone, bin_x, stim_LR, stim_prob_LR, init_Q, use_para)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ave_likeli_mean,~] = ...
    HMM_Hybrid_update_211031_max_all(para, Correct_trial, reward_trial, Chosen_trial, binary_tone, bin_x, stim_LR, stim_prob_LR, init_Q, use_para);
ave_likeli_mean
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function post_prob = Qlearning_ori(Q, choice, reward, alpha)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Control alpha
temp_Q = Q(:,choice) + alpha .* (reward - Q(:,choice)); 

if choice == 1
    post_prob = [temp_Q, (1-alpha) * Q(:,2)];
else
    post_prob = [(1-alpha) * Q(:,1), temp_Q];
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function post_prob = Transition_learning_ori(Q, reward, alpha)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Control alpha

post_prob = Q + alpha .* (reward - Q);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stim_prob = softmax(temp_stim,beta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear stime_prob
stim_prob(1) = exp(beta * temp_stim(1)) / (exp(beta * temp_stim(1)) + exp(beta * temp_stim(2)));
stim_prob(2) = exp(beta * temp_stim(2)) / (exp(beta * temp_stim(1)) + exp(beta * temp_stim(2)));

return


