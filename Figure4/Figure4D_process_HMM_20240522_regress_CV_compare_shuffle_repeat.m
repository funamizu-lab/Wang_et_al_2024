
%{
----------------------------------------------------------------------------
Determine the time window for task relevant neurons
%Start of task
%Sound on
%Sound off
%Before choice
%After choice (0sec)
%After choice (1sec)
%After choice (2sec)
%p = 0.001
%Each epoch, predict the prior, sensory and choice (integration)
----------------------------------------------------------------------------
%}
function Figure4D_process_HMM_20240522_regress_CV_compare_shuffle_repeat

[OFC_p, OFC_neuron] = temp_HMM_20240522_regress_CV_compare_shuffle('repeat_OFC_20230427');
[Hippo_p, Hippo_neuron] = temp_HMM_20240522_regress_CV_compare_shuffle('repeat_Hippo_20230427');
[AC_p, AC_neuron] = temp_HMM_20240522_regress_CV_compare_shuffle('repeat_AC_20230427');
[PPC_p, PPC_neuron] = temp_HMM_20240522_regress_CV_compare_shuffle('repeat_PPC_20230427');
close all

% length_neuron = [sig_before_sound, N_before_sound;
%                 sig_all_sound,     N_all_sound;
%                 sig_choice2,       N_choice2];
disp('plot before_sound')
plot_all_repeat(1, OFC_p,OFC_neuron,Hippo_p,Hippo_neuron,AC_p,AC_neuron,PPC_p,PPC_neuron);
disp('plot sound')
plot_all_repeat(2, OFC_p,OFC_neuron,Hippo_p,Hippo_neuron,AC_p,AC_neuron,PPC_p,PPC_neuron);
disp('plot choice')
plot_all_repeat(3, OFC_p,OFC_neuron,Hippo_p,Hippo_neuron,AC_p,AC_neuron,PPC_p,PPC_neuron);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_all_repeat(number, OFC_p,OFC_neuron,Hippo_p,Hippo_neuron,AC_p,AC_neuron,PPC_p,PPC_neuron)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y_ofc = [1/OFC_neuron(number,2):1/OFC_neuron(number,2):1]';
y_Hippo = [1/Hippo_neuron(number,2):1/Hippo_neuron(number,2):1]';
y_AC = [1/AC_neuron(number,2):1/AC_neuron(number,2):1]';
y_PPC = [1/PPC_neuron(number,2):1/PPC_neuron(number,2):1]';
% y_ofc = [1:-1/OFC_neuron(number,2):1/OFC_neuron(number,2)]';
% y_Hippo = [1:-1/Hippo_neuron(number,2):1/Hippo_neuron(number,2)]';
% y_AC = [1:-1/AC_neuron(number,2):1/AC_neuron(number,2)]';
% y_PPC = [1:-1/PPC_neuron(number,2):1/PPC_neuron(number,2)]';
% y_ofc = flipud(y_ofc);
% y_Hippo = flipud(y_Hippo);
% y_AC = flipud(y_AC);
% y_PPC = flipud(y_PPC);

OFC = sort(OFC_p(number).matrix); 
Hippo = sort(Hippo_p(number).matrix); 
AC = sort(AC_p(number).matrix); 
PPC = sort(PPC_p(number).matrix); 

OFC = remove_inf_OFC(OFC);
Hippo = remove_inf_OFC(Hippo);
AC = remove_inf_OFC(AC);
PPC = remove_inf_OFC(PPC);

thre = -log10(0.025);

figure
plot(OFC,y_ofc,'r')
hold on
plot(Hippo,y_Hippo,'g') 
hold on
plot(AC,y_AC,'b')
hold on
plot(PPC,y_PPC,'c')
hold on
plot([thre,thre],[0,1],'k')

[OFC_neuron(number,1:2); Hippo_neuron(number,1:2); AC_neuron(number,1:2); PPC_neuron(number,1:2)]

[OFC_neuron(number,3); Hippo_neuron(number,3); AC_neuron(number,3); PPC_neuron(number,3)]

return

% %before sound
% figure
% subplot(2,2,1)
% plot(OFC,y_ofc,'b')
% hold on
% plot([thre,thre],[0,1],'k')
% 
% subplot(2,2,2)
% plot(Hippo,y_Hippo,'b')
% hold on
% plot([thre,thre],[0,1],'k')
% 
% subplot(2,2,3)
% plot(AC,y_AC,'b')
% hold on
% plot([thre,thre],[0,1],'k')
% 
% subplot(2,2,4)
% plot(PPC,y_PPC,'b')
% hold on
% plot([thre,thre],[0,1],'k')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OFC = remove_inf_OFC(OFC)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp = find(OFC == inf);
temp_OFC = OFC;
temp_OFC(temp) = [];
max_OFC = max(temp_OFC);
OFC(temp) = max_OFC;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p_values, length_neuron] = temp_HMM_20240522_regress_CV_compare_shuffle(folders)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch nargin
    case 0
        hoge
    case 1
        disp('OK to analyze')
    otherwise
        hoge
end

[analysis_dir,depth_def] = eval(folders);

prob_before_sound = [];
prob_all_sound = [];
prob_choice2 = [];
for i = 1:length(analysis_dir)
    [i,length(analysis_dir)]
    temp_dir = analysis_dir{i};
    
    [temp_before_sound,temp_all_sound,temp_choice] = ...
        HMM_ephys_20240522_regress_get_data_shuffle(temp_dir, depth_def);

    prob_before_sound = [prob_before_sound; temp_before_sound];
    prob_all_sound = [prob_all_sound; temp_all_sound];
    prob_choice2 = [prob_choice2; temp_choice];
end
delete(gcp('nocreate'))

%check whether the STD = 1.96 is about 0.025
temp = find(prob_before_sound(:,1) > 0.01)
prob_before_sound(temp,:)
figure
plot(prob_before_sound(temp,1),prob_before_sound(temp,2),'k.')

figure
subplot(1,3,1)
plot(prob_before_sound(:,1), prob_before_sound(:,2), 'b.')
subplot(1,3,2)
plot(prob_all_sound(:,1), prob_all_sound(:,2), 'b.')
subplot(1,3,3)
plot(prob_choice2(:,1), prob_choice2(:,2), 'b.')

figure
subplot(1,3,1)
plot(-log10(prob_before_sound(:,1)), prob_before_sound(:,2), 'r.')
subplot(1,3,2)
plot(-log10(prob_all_sound(:,1)), prob_all_sound(:,2), 'r.')
subplot(1,3,3)
plot(-log10(prob_choice2(:,1)), prob_choice2(:,2), 'r.')
%convert to p_values
% sum(prob_before_sound)
% sum(prob_all_sound)
% sum(prob_choice2)

[N_before_sound,~] = size(prob_before_sound)
[N_all_sound,~] = size(prob_all_sound)
[N_choice2,~] = size(prob_choice2)

%We need to plot the -log10(:,1)
p_values(1).matrix = -log10(prob_before_sound(:,1));
p_values(2).matrix = -log10(prob_all_sound(:,1));
p_values(3).matrix = -log10(prob_choice2(:,1));

thre = 0.025;
% sig_before_sound = length(find(prob_before_sound(:,1) < thre));
% sig_all_sound = length(find(prob_all_sound(:,1) < thre));
% sig_choice2 = length(find(prob_choice2(:,1) < thre));
sig_before_sound = length(find(p_values(1).matrix > -log10(thre)));
sig_all_sound = length(find(p_values(2).matrix > -log10(thre)));
sig_choice2 = length(find(p_values(3).matrix > -log10(thre)));

length_neuron = [sig_before_sound, N_before_sound, sig_before_sound/N_before_sound;
                sig_all_sound,     N_all_sound,    sig_all_sound/N_all_sound;
                sig_choice2,       N_choice2,      sig_choice2/N_choice2];
length_neuron

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [use_prob_b, use_prob_p, size_neuron] = analyze_b_regress(use_b, use_p)
%1:constant
%2:choice
%3:sound
%4:correct/error
%5: previous_choice
%6: previous_sound
%7: previous_correct_error
%8: mouse velocity
[size_regress,size_neuron] = size(use_b);

for i = 1:size_regress
    temp_b = use_b(i,:);
    temp_p = use_p(i,:);
    
    temp_b = length(find(temp_b ~= 0));
    temp_p = length(find(temp_p < 0.01));
    prob_b(i) = temp_b ./ size_neuron;
    prob_p(i) = temp_p ./ size_neuron;
end

use_prob_b = prob_b(2:8);
use_prob_p = prob_p(2:8);

return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [PostLow_left, PostHigh_left, PostLow_right, PostHigh_right] = ...
    test_surprise_choice(all_norm_spike, low_sig_neuron)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PostLow_left = all_norm_spike(low_sig_neuron,15);
PostHigh_left = all_norm_spike(low_sig_neuron,16);
PostLow_right = all_norm_spike(low_sig_neuron,17);
PostHigh_right = all_norm_spike(low_sig_neuron,18);

figure
subplot(1,2,1)
plot(PostLow_left, PostLow_right, 'b.')
hold on
plot([-1, 10],[-1, 10],'k')
set(gca,'xlim',[-1, 10],'ylim',[-1, 10])
subplot(1,2,2)
plot(PostHigh_left, PostHigh_right, 'r.')
hold on
plot([-1, 10],[-1, 10],'k')
set(gca,'xlim',[-1, 10],'ylim',[-1, 10])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function plot_surprise_choice(PostLow_left, PostHigh_left, PostLow_right, PostHigh_right)
function plot_surprise_choice(Repeat_prefer, Switch_prefer, Switch_nonprefer, Repeat_nonprefer)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
predict_low = Repeat_prefer - Switch_prefer;
predict_high = Switch_nonprefer - Repeat_nonprefer;
choice_low = Repeat_prefer - Switch_nonprefer;
choice_high = Switch_prefer - Repeat_nonprefer;
predict_ave = mean([predict_low, predict_high],2);
choice_ave = mean([choice_low, choice_high],2);

% figure
% plot(predict_ave, choice_ave, 'k.')

figure
subplot(2,3,1)
plot(Repeat_prefer, Switch_prefer, 'b.') %same choice, X_repeat or Y_zigzag
hold on
plot([-1, 10],[-1, 10],'k')
set(gca,'xlim',[-1, 10],'ylim',[-1, 10])
subplot(2,3,2)
plot(Switch_nonprefer, Repeat_nonprefer, 'r.')  %same choice, X_zigzag or Y_repeat
hold on
plot([-1, 10],[-1, 10],'k')
set(gca,'xlim',[-1, 10],'ylim',[-1, 10])

subplot(2,3,3)
boxplot([Repeat_prefer-Switch_prefer, Switch_nonprefer-Repeat_nonprefer])

subplot(2,3,4)
plot(Repeat_prefer, Switch_nonprefer, 'b.')
hold on
plot([-1, 10],[-1, 10],'k')
set(gca,'xlim',[-1, 10],'ylim',[-1, 10])
subplot(2,3,5)
plot(Switch_prefer, Repeat_nonprefer, 'r.')
hold on
plot([-1, 10],[-1, 10],'k')
set(gca,'xlim',[-1, 10],'ylim',[-1, 10])

subplot(2,3,6)
boxplot([Repeat_prefer-Switch_nonprefer, Switch_prefer-Repeat_nonprefer])

[median(Repeat_prefer), median(Switch_prefer)]
signrank(Repeat_prefer, Switch_prefer)

[median(Switch_nonprefer), median(Repeat_nonprefer)]
signrank(Switch_nonprefer, Repeat_nonprefer)

[length(Repeat_prefer),length(Switch_prefer),length(Switch_nonprefer), length(Repeat_nonprefer)]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [PostLow_low, PostHigh_low, PostLow_high, PostHigh_high] = ...
%     test_surprise_sound(all_norm_spike, low_sig_neuron)
function test_surprise_sound(all_norm_spike, low_sig_neuron)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PostLow = previous tone was low, _low is current is low
PostLow_low = all_norm_spike(low_sig_neuron,19);
PostHigh_low = all_norm_spike(low_sig_neuron,20);
PostLow_high = all_norm_spike(low_sig_neuron,21);
PostHigh_high = all_norm_spike(low_sig_neuron,22);

figure
subplot(2,2,1)
plot(PostLow_low, PostHigh_low, 'b.') 
hold on
plot([-1, 10],[-1, 10],'k')
set(gca,'xlim',[-1, 10],'ylim',[-1, 10])
subplot(2,2,2)
plot(PostLow_high, PostHigh_high, 'r.')
hold on
plot([-1, 10],[-1, 10],'k')
set(gca,'xlim',[-1, 10],'ylim',[-1, 10])

subplot(2,2,3)
plot(PostLow_low, PostLow_high, 'b.')
hold on
plot([-1, 10],[-1, 10],'k')
set(gca,'xlim',[-1, 10],'ylim',[-1, 10])
subplot(2,2,4)
plot(PostHigh_low, PostHigh_high, 'r.')
hold on
plot([-1, 10],[-1, 10],'k')
set(gca,'xlim',[-1, 10],'ylim',[-1, 10])

[median(PostLow_low), median(PostHigh_low)]
signrank(PostLow_low, PostHigh_low)

[median(PostLow_high), median(PostHigh_high)]
signrank(PostLow_high, PostHigh_high)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function neuron_sig = hist_index_sound(index1, p1, thre)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sig_neuron = find(p1 < thre);
sig_index1 = index1(sig_neuron);

neuron_sig = [length(sig_index1), length(index1)]

hist_x = [-1:0.1:1];
plot_hist_x = hist_x(1:end-1) + 0.05;

hist_all = histcounts(index1, hist_x);
hist_sig = histcounts(sig_index1, hist_x);

plot(plot_hist_x, hist_all,'k')
hold on
area(plot_hist_x, hist_sig, 'FaceColor',[0.5 0.5 0.5],'EdgeColor','none')
set(gca,'xlim',[-1 1])

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p, median_value,rho,pval] = plot_index_abs_index_sig(index1, index2, right_neuron,left_neuron,non_sig_neuron)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
subplot(1,2,1)
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

[rho,pval] = corr(index1,index2);

subplot(1,2,2)
plot([0 1],[0 1],'k:')
hold on
%plot(abs(index1), abs(index2), 'b.')
plot(abs(index1(non_sig_neuron)), abs(index2(non_sig_neuron)), 'k.')
hold on
plot(abs(index1(left_neuron)), abs(index2(left_neuron)), 'b.')
hold on
plot(abs(index1(right_neuron)), abs(index2(right_neuron)), 'r.')

p = signrank(abs(index1), abs(index2));
median_value = [median(abs(index1)), median(abs(index2))];

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p, median_value,rho,pval] = plot_index_abs_index(index1, index2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
subplot(1,2,1)
plot([-1 1],[0 0],'k:')
hold on
plot([0 0],[-1 1],'k:')
hold on
plot([-1 1],[-1 1],'k:')
hold on
plot(index1, index2, 'b.')

[rho,pval] = corr(index1,index2);

subplot(1,2,2)
plot([0 1],[0 1],'k:')
hold on
plot(abs(index1), abs(index2), 'b.')

p = signrank(abs(index1), abs(index2));
median_value = [median(abs(index1)), median(abs(index2))];

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_basic_tone_choice_index(Sound_prefer_index, use_neuron)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
subplot(1,3,1) %Short sound
plot(Sound_prefer_index(use_neuron,4), Sound_prefer_index(use_neuron,1),'b.')
hold on
plot([0 0],[-1 1],'k--')
hold on
plot([-1 1],[0 0],'k--')
set(gca,'xlim',[-1 1],'ylim',[-1 1])
%signrank(Sound_prefer_index(:,4), Sound_prefer_index(:,1))
subplot(1,3,2) %Long sound
plot(Sound_prefer_index(use_neuron,5), Sound_prefer_index(use_neuron,2),'b.')
hold on
plot([0 0],[-1 1],'k--')
hold on
plot([-1 1],[0 0],'k--')
set(gca,'xlim',[-1 1],'ylim',[-1 1])
%signrank(Sound_prefer_index(:,5), Sound_prefer_index(:,2))
subplot(1,3,3) %Long sound end
plot(Sound_prefer_index(use_neuron,6), Sound_prefer_index(use_neuron,3),'b.')
hold on
plot([0 0],[-1 1],'k--')
hold on
plot([-1 1],[0 0],'k--')
set(gca,'xlim',[-1 1],'ylim',[-1 1])
%signrank(Sound_prefer_index(:,6), Sound_prefer_index(:,3))

signrank(Sound_prefer_index(:,4))
signrank(Sound_prefer_index(:,5))
signrank(Sound_prefer_index(:,6))

figure %Onset and Offset of sound @ correct trials
plot(Sound_prefer_index(use_neuron,3),Sound_prefer_index(use_neuron,2),'b.')
hold on
plot([-1 1],[-1 1],'k')
hold on
plot([0 0],[-1 1],'k--')
hold on
plot([-1 1],[0 0],'k--')
set(gca,'xlim',[-1 1],'ylim',[-1 1])
signrank(Sound_prefer_index(use_neuron,2), Sound_prefer_index(use_neuron,3))
signrank(abs(Sound_prefer_index(use_neuron,2)), abs(Sound_prefer_index(use_neuron,3)))
[median(Sound_prefer_index(use_neuron,2)), median(Sound_prefer_index(use_neuron,3))]

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_sigR_sig_L(sig_all,sig_R,sig_L,tuning0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sig_all1=sig_all(1).matrix(tuning0,:);
sig_all2=sig_all(2).matrix(tuning0,:);
sig_all3=sig_all(3).matrix(tuning0,:);
sig_R1=sig_R(1).matrix(tuning0,:);
sig_R2=sig_R(2).matrix(tuning0,:);
sig_R3=sig_R(3).matrix(tuning0,:);
sig_L1=sig_L(1).matrix(tuning0,:);
sig_L2=sig_L(2).matrix(tuning0,:);
sig_L3=sig_L(3).matrix(tuning0,:);

figure
subplot(1,3,1)
plot(nanmean(sig_all1),'k')
hold on
plot(nanmean(sig_R1),'r')
hold on
plot(nanmean(sig_L1),'b')
subplot(1,3,2)
plot(nanmean(sig_all2),'k')
hold on
plot(nanmean(sig_R2),'r')
hold on
plot(nanmean(sig_L2),'b')
subplot(1,3,3)
plot(nanmean(sig_all3),'k')
hold on
plot(nanmean(sig_R3),'r')
hold on
plot(nanmean(sig_L3),'b')

clear p
for i = 1:6
    p(1,i) = signrank(sig_R1(:,i),sig_L1(:,i));
    p(2,i) = signrank(sig_R2(:,i),sig_L2(:,i));
    p(3,i) = signrank(sig_R3(:,i),sig_L3(:,i));
end
p