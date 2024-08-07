
function FigureS5E_CurrentOutcome_20240703_surprise_before_sound_depth

% % % % % repeat
[repeat_OFC_prefer_correct, repeat_OFC_nonprefer_correct, repeat_OFC_prefer_error, repeat_OFC_nonprefer_error, ...
    repeat_p(1,:), repeat_neuron(1)] = S5Eprocess_CurrentOutcome_20240703_surprise_before_sound_depth('repeat_OFC_20240427');

[repeat_AC_prefer_correct, repeat_AC_nonprefer_correct, repeat_AC_prefer_error, repeat_AC_nonprefer_error, ...
    repeat_p(2,:), repeat_neuron(2)] = S5Eprocess_CurrentOutcome_20240703_surprise_before_sound_depth('repeat_AC_20240427');

% % % % % % zigzag
[zigzag_OFC_prefer_correct, zigzag_OFC_nonprefer_correct, zigzag_OFC_prefer_error, zigzag_OFC_nonprefer_error,...
    zigzag_p(1,:), zigzag_neuron(1)] = S5Eprocess_CurrentOutcome_20240703_surprise_before_sound_depth('altern_OFC_20240427');

[zigzag_AC_prefer_correct, zigzag_AC_nonprefer_correct, zigzag_AC_prefer_error, zigzag_AC_nonprefer_error,...
    zigzag_p(2,:), zigzag_neuron(2)] = S5Eprocess_CurrentOutcome_20240703_surprise_before_sound_depth('altern_AC_20240427');

[zigzag_PPC_prefer_correct, zigzag_PPC_nonprefer_correct, zigzag_PPC_prefer_error, zigzag_PPC_nonprefer_error,...
    zigzag_p(3,:), zigzag_neuron(3)] = S5Eprocess_CurrentOutcome_20240703_surprise_before_sound_depth('altern_PPC_20240427');

[zigzag_M1_prefer_correct, zigzag_M1_nonprefer_correct, zigzag_M1_prefer_error, zigzag_M1_nonprefer_error,...
    zigzag_p(4,:), zigzag_neuron(4)] = S5Eprocess_CurrentOutcome_20240703_surprise_before_sound_depth('altern_M1_20240427');

[zigzag_STR_prefer_correct, zigzag_STR_nonprefer_correct, zigzag_STR_prefer_error, zigzag_STR_nonprefer_error,...
    zigzag_p(5,:), zigzag_neuron(5)] = S5Eprocess_CurrentOutcome_20240703_surprise_before_sound_depth('altern_STR_20240427');

close all


repeat_prefer_OFC =[repeat_OFC_prefer_correct,repeat_OFC_prefer_error];
repeat_nonprefer_OFC =[repeat_OFC_nonprefer_correct,repeat_OFC_nonprefer_error];
repeat_prefer_AC =[repeat_AC_prefer_correct,repeat_AC_prefer_error];
repeat_nonprefer_AC =[repeat_AC_nonprefer_correct,repeat_AC_nonprefer_error];

zigzag_prefer_OFC =[zigzag_OFC_prefer_correct,zigzag_OFC_prefer_error];
zigzag_nonprefer_OFC =[zigzag_OFC_nonprefer_correct,zigzag_OFC_nonprefer_error];
zigzag_prefer_AC =[zigzag_AC_prefer_correct,zigzag_AC_prefer_error];
zigzag_nonprefer_AC =[zigzag_AC_nonprefer_correct,zigzag_AC_nonprefer_error];
zigzag_prefer_PPC =[zigzag_PPC_prefer_correct,zigzag_PPC_prefer_error];
zigzag_nonprefer_PPC =[zigzag_PPC_nonprefer_correct,zigzag_PPC_nonprefer_error];
zigzag_prefer_M1 =[zigzag_M1_prefer_correct,zigzag_M1_prefer_error];
zigzag_nonprefer_M1 =[zigzag_M1_nonprefer_correct,zigzag_M1_nonprefer_error];
zigzag_prefer_STR =[zigzag_STR_prefer_correct,zigzag_STR_prefer_error];
zigzag_nonprefer_STR =[zigzag_STR_nonprefer_correct,zigzag_STR_nonprefer_error];


temp1 = cell(2,4);
temp2 = cell(2,10);

AArep = cell(1,2);BBrep = cell(1,2);CCrep = cell(1,2);DDrep = cell(1,2);
for ii=1:size(temp1,1)
    AArep{ii}=repeat_prefer_OFC(:,ii);
    BBrep{ii}=repeat_prefer_AC(:,ii);
    CCrep{ii}=repeat_nonprefer_OFC(:,ii);
    DDrep{ii}=repeat_nonprefer_AC(:,ii);
end
AAzig = cell(1,2);BBzig = cell(1,2);CCzig = cell(1,2);DDzig = cell(1,2);EEzig = cell(1,2);
FFzig = cell(1,2);GGzig = cell(1,2);HHzig = cell(1,2);IIzig = cell(1,2);JJzig = cell(1,2);
for ii=1:size(temp2,1)
    AAzig{ii}=zigzag_prefer_OFC(:,ii);
    BBzig{ii}=zigzag_prefer_AC(:,ii);
    CCzig{ii}=zigzag_prefer_PPC(:,ii);
    DDzig{ii}=zigzag_prefer_M1(:,ii);
    EEzig{ii}=zigzag_prefer_STR(:,ii);
    FFzig{ii}=zigzag_nonprefer_OFC(:,ii);
    GGzig{ii}=zigzag_nonprefer_AC(:,ii);
    HHzig{ii}=zigzag_nonprefer_PPC(:,ii);
    IIzig{ii}=zigzag_nonprefer_M1(:,ii);
    JJzig{ii}=zigzag_nonprefer_STR(:,ii);
end

data_repeat=vertcat(AArep,BBrep,CCrep,DDrep);
data_zigzag = vertcat(AAzig,BBzig,CCzig,DDzig,EEzig,FFzig,GGzig,HHzig,IIzig,JJzig);

Mlab = {'Current correct', 'Current error'};
colorcorrect = [0.4660 0.6740 0.1880]';
colorerror = [0.6350 0.0780 0.1840]';
color = [colorerror,colorcorrect];

figure
multiple_boxplot(data_repeat,{'OFC','AC','OFC','AC'},Mlab,color)
hold on
yline(0,':k')
ylim([-4 4])
legend(fliplr(Mlab),'Location','southeast');
title('Repeating condition (Before sound)')
text(1.5, 3.3, 'Preferred Choice')
text(2.8, 3.3, 'Non-preferred Choice')
set(gcf,'Position',[441,351,470,360])


figure
multiple_boxplot(data_zigzag,{'OFC','AC','PPC','M1','STR','OFC','AC','PPC','M1','STR'},Mlab,color)
hold on
yline(0,':k')
ylim([-4 4])
legend(fliplr(Mlab),'Location','southeast');
title('Alternating condition (Before sound)')
text(2, 3, 'Preferred Choice')
text(6.5, 3, 'Non-preferred Choice')
set(gcf,'Position',[389,433,920,360])

return


function multiple_boxplot(data,xlab,~,colors)


% check that data is ok.
if ~iscell(data)
    error('Input data is not even a cell array!');
end

% Get sizes
M=size(data,2);
L=size(data,1);
if nargin>=4
    if size(colors,2)~=M
        error('Wrong amount of colors!');
    end
end
if nargin>=2
    if length(xlab)~=L
        error('Wrong amount of X labels given');
    end
end

% Calculate the positions of the boxes
positions=1:0.25:M*L*0.25+1+0.25*L;
positions(1:M+1:end)=[];

% Extract data and label it in the group correctly
x=[];
group=[];
for ii=1:L
    for jj=1:M
        aux=data{ii,jj};
        x=vertcat(x,aux(:));
        group=vertcat(group,ones(size(aux(:)))*jj+(ii-1)*M);
    end
end
% Plot it

boxplot(x,group, 'positions', positions,'Symbol',' ','BoxStyle','outline');

% Set the Xlabels
aux=reshape(positions,M,[]);
labelpos = sum(aux,1)./M;

set(gca,'xtick',labelpos)
if nargin>=2
    set(gca,'xticklabel',xlab);
else
    idx=1:L;
    set(gca,'xticklabel',strsplit(num2str(idx),' '));
end
    

% Get some colors
if nargin>=4
    cmap=colors;
else
    cmap = hsv(M);
    cmap=vertcat(cmap',ones(1,M)*0.5);
end
color=repmat(cmap, 1, L);

% Apply colors
h = findobj(gca,'Tag','Box');
for jj=1:length(h)
   patch(get(h(jj),'XData'),get(h(jj),'YData'),color(1:3,jj)','EdgeColor',color(1:3,jj),'FaceColor','none');
end

% if nargin>=3
%     legend(fliplr(Mlab));
% end
return


function [prefer_sabun_correct, nonprefer_sabun_correct, prefer_sabun_error, nonprefer_sabun_error,...
    p_sabun, length_neuron] ...
    = S5Eprocess_CurrentOutcome_20240703_surprise_before_sound_depth(folders, kaiseki_number)

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

Brainarea = [folders(1:6),' ', folders(8:10)];
[analysis_dir,depth_def] = eval(folders);

all_pre_sound = [];
all_p_pre_sound = [];
all_norm_spike = [];
all_p_surprise = [];
all_pre_choice_correct_error = [];
all_p_pre_choice_correct_error = [];

for i = 1:length(analysis_dir)
    disp([i,length(analysis_dir)])
    temp_dir = analysis_dir{i};
    
    [pre_sound, p_pre_sound, norm_spike, ~, ...
    ~, ~, ~,...
    ~, ~, ~, ~,...
    pre_choice_correct_error, p_pre_choice_correct_error,~,~,~,~,~,~,~,~,~,~,~,~,~,~] = ...
        S5E_HMM_ephys_20240401_surprise_before_sound_depth(temp_dir, kaiseki_number,depth_def);

    all_pre_sound = [all_pre_sound; pre_sound];
    all_p_pre_sound = [all_p_pre_sound; p_pre_sound];
    all_norm_spike = [all_norm_spike; norm_spike];

    all_pre_choice_correct_error = [all_pre_choice_correct_error; pre_choice_correct_error]; % pre choice correct index
    all_p_pre_choice_correct_error = [all_p_pre_choice_correct_error; p_pre_choice_correct_error]; % pre choice error index
end
delete(gcp('nocreate'))

disp([length(all_pre_sound),length(all_norm_spike),length(all_p_surprise)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot about low, high neurons
% right_neuron = find(all_neuron_choice_category == 1);
% left_neuron = find(all_neuron_choice_category == 0);

if contains(Brainarea,'repeat')
    right_neuron = find(all_pre_choice_correct_error(:,1) >= 0);
    left_neuron = find(all_pre_choice_correct_error(:,1) < 0);
    sig_sound = find(all_p_pre_sound(:,3) < 0.01); %Sig diff. between previous low and high
else
    left_neuron = find(all_pre_choice_correct_error(:,1) >= 0);
    right_neuron = find(all_pre_choice_correct_error(:,1) < 0);
    sig_sound = find(all_p_pre_sound(:,3) < 0.01); %Sig diff. between previous low and high

end

right_sig_neuron = intersect(right_neuron, sig_sound);
left_sig_neuron = intersect(left_neuron, sig_sound);
non_sig_neuron = setdiff(1:length(all_p_pre_sound), sig_sound);
check_neuron = length(right_sig_neuron) + length(left_sig_neuron) + length(non_sig_neuron);
if check_neuron ~= length(all_p_pre_sound)
    hoge
else
    disp([length(left_sig_neuron),length(right_sig_neuron),length(non_sig_neuron)])
end

%Use only on low sig neuron
[L_sig_PostLow_leftcorrect, L_sig_PostHigh_leftcorrect, L_sig_PostLow_rightcorrect, L_sig_PostHigh_rightcorrect,...
    L_sig_PostLow_lefterror, L_sig_PostHigh_lefterror, L_sig_PostLow_righterror, L_sig_PostHigh_righterror] = ...
    test_surprise_choice(all_norm_spike, left_sig_neuron);
%Use only on high sig neuron
[R_sig_PostLow_leftcorrect, R_sig_PostHigh_leftcorrect, R_sig_PostLow_rightcorrect, R_sig_PostHigh_rightcorrect,....
    R_sig_PostLow_lefterror, R_sig_PostHigh_lefterror, R_sig_PostLow_righterror, R_sig_PostHigh_righterror] = ...
    test_surprise_choice(all_norm_spike, right_sig_neuron);

%Activity is already flipped

if contains(Brainarea,'repeat')
    Repeat_prefer_correct = [L_sig_PostLow_leftcorrect; R_sig_PostHigh_rightcorrect];
    Repeat_nonprefer_correct = [L_sig_PostHigh_rightcorrect; R_sig_PostLow_leftcorrect];
    Switch_prefer_correct = [L_sig_PostLow_rightcorrect; R_sig_PostHigh_leftcorrect];
    Switch_nonprefer_correct = [L_sig_PostHigh_leftcorrect; R_sig_PostLow_rightcorrect];

    Repeat_prefer_error = [L_sig_PostLow_lefterror; R_sig_PostHigh_righterror];
    Repeat_nonprefer_error = [L_sig_PostHigh_righterror; R_sig_PostLow_lefterror];
    Switch_prefer_error = [L_sig_PostLow_righterror; R_sig_PostHigh_lefterror];
    Switch_nonprefer_error = [L_sig_PostHigh_lefterror; R_sig_PostLow_righterror];
else
    
    Repeat_prefer_correct = [L_sig_PostHigh_rightcorrect; R_sig_PostLow_leftcorrect];
    Repeat_nonprefer_correct = [L_sig_PostLow_leftcorrect; R_sig_PostHigh_rightcorrect];
    Switch_prefer_correct = [L_sig_PostHigh_leftcorrect; R_sig_PostLow_rightcorrect];
    Switch_nonprefer_correct = [L_sig_PostLow_rightcorrect; R_sig_PostHigh_leftcorrect];


    Repeat_prefer_error = [L_sig_PostHigh_righterror; R_sig_PostLow_lefterror];
    Repeat_nonprefer_error = [L_sig_PostLow_lefterror; R_sig_PostHigh_righterror];
    Switch_prefer_error = [L_sig_PostHigh_lefterror; R_sig_PostLow_righterror];
    Switch_nonprefer_error = [L_sig_PostLow_righterror; R_sig_PostHigh_lefterror];


   
end

p_prefer_correct = signrank(Repeat_prefer_correct,Switch_prefer_correct);
p_nonprefer_correct = signrank(Repeat_nonprefer_correct,Switch_nonprefer_correct);


p_prefer_error = signrank(Repeat_prefer_error,Switch_prefer_error);
p_nonprefer_error = signrank(Repeat_nonprefer_error,Switch_nonprefer_error);

prefer_sabun_correct = Repeat_prefer_correct - Switch_prefer_correct;
prefer_sabun_error = Repeat_prefer_error - Switch_prefer_error;

nonprefer_sabun_correct = Repeat_nonprefer_correct - Switch_nonprefer_correct;
nonprefer_sabun_error = Repeat_nonprefer_error - Switch_nonprefer_error;

p_sabun = [p_prefer_correct, p_nonprefer_correct, p_prefer_error, p_nonprefer_error];

temp1 = length(prefer_sabun_correct);
temp2 = length(nonprefer_sabun_correct);

if temp1 ~= temp2
    hoge
else
    length_neuron = temp1;
end


return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [PostLow_Leftcorrect, PostHigh_Leftcorrect, PostLow_Rightcorrect, PostHigh_Rightcorrect,...
    PostLow_Lefterror,PostHigh_Lefterror,PostLow_Righterror,PostHigh_Righterror] = ...
    test_surprise_choice(all_norm_spike, low_sig_neuron)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PostLow_Leftcorrect = all_norm_spike(low_sig_neuron,33); %PostLowCorrect_left
PostHigh_Leftcorrect = all_norm_spike(low_sig_neuron,37); %PostHighCorrect_left
PostLow_Rightcorrect = all_norm_spike(low_sig_neuron,35); %PostLowCorrect_right
PostHigh_Rightcorrect = all_norm_spike(low_sig_neuron,39); %PostHighCorrect_right

PostLow_Lefterror = all_norm_spike(low_sig_neuron,34); %PostLowCorrect_left
PostHigh_Lefterror = all_norm_spike(low_sig_neuron,38); %PostHighCorrect_left
PostLow_Righterror = all_norm_spike(low_sig_neuron,36); %PostLowCorrect_right
PostHigh_Righterror = all_norm_spike(low_sig_neuron,40); %PostHighCorrect_right

function plot_surprise_choice_prefer(Repeat_prefer, Switch_prefer, plot_color)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% figure
% plot(predict_ave, choice_ave, 'k.')
plot(Repeat_prefer, Switch_prefer, '.','color', plot_color) %same choice, X_repeat or Y_zigzag
hold on
plot([-1, 10],[-1, 10],'k')
set(gca,'xlim',[-1, 10],'ylim',[-1, 10])


function plot_surprise_choice_nonprefer(Switch_nonprefer, Repeat_nonprefer, plot_color)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% figure
% plot(predict_ave, choice_ave, 'k.')

plot(Switch_nonprefer, Repeat_nonprefer, '.','color', plot_color)  %same choice, X_zigzag or Y_repeat
hold on
plot([-1, 10],[-1, 10],'k')
set(gca,'xlim',[-1, 10],'ylim',[-1, 10])


%{
----------------------------------------------------------------------------
First_take number of frames in each tif files
Analyzing imaging data simply
At least for the correct rate
----------------------------------------------------------------------------
%}
function [pre_sound, p_pre_sound, norm_spike, p_surprise, ...
    choice_correct_error, p_choice_correct_error, neuron_choice_category,...
    tone_correct_error, p_tone_correct_error, pre_tone_correct_error, p_pre_tone_correct_error,...
    pre_choice_correct_error, p_pre_choice_correct_error,...
    SameDif_preTone_correct_error,p_SameDif_preTone_correct_error,...
    preCorrectTone_index_correct_correct_error,p_preCorrectTone_index_correct_correct_error,...
    preCorrectChoice_index_correct_correct_error,p_preCorrectChoice_index_correct_correct_error,...
    pre_same_pLcL_pLcH_index,p_pre_same_pLcL_pLcH_index,...
    pre_same_pHcL_pHcH_index,p_pre_same_pHcL_pHcH_index,...
    cur_same_pLcL_pHcL_index,p_cur_same_pLcL_pHcL_index,...
    cur_same_pLcH_pHcH_index,p_cur_same_pLcH_pHcH_index] = ...
    S5E_HMM_ephys_20240401_surprise_before_sound_depth(pathname, kaiseki_number,depth_def)


switch nargin
    case 0
        pathname = pwd;
    case 3
        disp('OK to analyze')
    otherwise
        hoge
end
cd(pathname)

% [filename1, pathname1,findex]=uigetfile('*.*','frame file');
% filename1 = [pathname1,filename1];
% load(filename1)
temp = dir('HMM_spike_count_neurons_sound_trans_20240326*');
if length(temp) ~= 1
    hoge
end
load(temp.name);
%neuron_index p_index

temp = dir('sig_HMM_neurons_20230310*');
if length(temp) ~= 1
    hoge
end
load(temp.name);
%p_task: around sound
%p_task2: around choice

temp = dir('depth_spike_20230427*');
if length(temp) == 1
    load(temp.name);
    %spike_depth def_depth length_neuron
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
    hoge
end

% %Get sig neurons
new_p_thre = 10;
% [size_neuron,~] = size(p_task);
% sig1 = get_sig_neuron_all(p_task,new_p_thre);
% sig2 = get_sig_neuron_all(p_task2,new_p_thre);
% sig_neuron = union(sig1,sig2);

%sig_before_sound = get_sig_neuron_time_window(p_task,new_p_thre,[6:15]); %1000ms before sound
sig_before_sound = get_sig_neuron_time_window(p_task,new_p_thre,10:15); %600ms before sound
sig_during_sound = get_sig_neuron_time_window(p_task,new_p_thre,16:21); %During sound
%sig_during_choice = get_sig_neuron_time_window(p_task2,new_p_thre,[6:8]); %During choice
sig_during_choice = get_sig_neuron_time_window(p_task2,new_p_thre,6:15); %During choice 1000ms

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%neuron_index.before_sound
%neuron_index.init_sound
%neuron_index.late_sound
%neuron_index.all_sound
%neuron_index.before_choice
%neuron_index.choice1
%neuron_index.choice2
%neuron_index.choice3
if kaiseki_number == 1 %Before sound
    use_sig_neuron = sig_before_sound;
    use_neuron_index = neuron_index.before_sound;
    use_p_index = p_index.before_sound;
    use_norm_spike = norm_spike_group.before_sound;
    use_p_surprise = p_surprise.before_sound;

elseif kaiseki_number == 2 %during sound
    use_sig_neuron = sig_during_sound;
    use_neuron_index = neuron_index.all_sound;
    use_p_index = p_index.all_sound;
    use_norm_spike = norm_spike_group.all_sound;
    use_p_surprise = p_surprise.all_sound;

elseif kaiseki_number == 3 %during choice %1000ms
    use_sig_neuron = sig_during_choice;
    use_neuron_index = neuron_index.choice2;
    use_p_index = p_index.choice2;
    use_norm_spike = norm_spike_group.choice2;
    use_p_surprise = p_surprise.choice2;

else
    hopge
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ADD depth definition
use_sig_neuron = intersect(use_sig_neuron,depth_neuron);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
%%%
%%%%%
%The activity was changed with previous choice and outcome, then how?
%Candidates:
%current choice
%previous low or high
%previous left or hight
%%%%%
%%%
%

%Choice index current
% 1. current choice index
% 2. current tone index
% 3. current outcome index
% 4. current choice correct index
% 5. current choice error index
% 6. pre tone index
% 7. pre choice index
% 12. tone correct index
% 13. tone error index
% 14. pre tone correct index
% 15. pre tone error index
% 16. pre choice correct index
% 17. pre choice error index


pre_sound(:,1) = use_neuron_index(use_sig_neuron,1); %current_choice
pre_sound(:,2) = use_neuron_index(use_sig_neuron,6); %pre sound
pre_sound(:,3) = use_neuron_index(use_sig_neuron,7); %pre choice
pre_sound(:,4) = use_neuron_index(use_sig_neuron,10); %pre SoundCorrect, current left
pre_sound(:,5) = use_neuron_index(use_sig_neuron,11); %pre SoundCorrect, current right

pre_sound(:,6) = use_neuron_index(use_sig_neuron,2); %current tone index


p_pre_sound(:,1) = use_p_index(use_sig_neuron,1); %current choice index
p_pre_sound(:,2) = use_p_index(use_sig_neuron,6); %pre tone index
p_pre_sound(:,3) = use_p_index(use_sig_neuron,7); %pre choice index
p_pre_sound(:,4) = use_p_index(use_sig_neuron,10); %pre SoundCorrect, current left
p_pre_sound(:,5) = use_p_index(use_sig_neuron,11); %pre SoundCorrect, current right

p_pre_sound(:,6) = use_p_index(use_sig_neuron,2); %current tone index

choice_correct_error(:,1) = use_neuron_index(use_sig_neuron,4); %Correct choice index
choice_correct_error(:,2) = use_neuron_index(use_sig_neuron,5); %Error choice index
p_choice_correct_error(:,1) = use_p_index(use_sig_neuron,4);
p_choice_correct_error(:,2) = use_p_index(use_sig_neuron,5);

tone_correct_error(:,1) = use_neuron_index(use_sig_neuron,12); %Correct choice index
tone_correct_error(:,2) = use_neuron_index(use_sig_neuron,13); %Error choice index
p_tone_correct_error(:,1) = use_p_index(use_sig_neuron,12);
p_tone_correct_error(:,2) = use_p_index(use_sig_neuron,13);

pre_tone_correct_error(:,1) = use_neuron_index(use_sig_neuron,14); %Correct choice index
pre_tone_correct_error(:,2) = use_neuron_index(use_sig_neuron,15); %Error choice index
p_pre_tone_correct_error(:,1) = use_p_index(use_sig_neuron,14);
p_pre_tone_correct_error(:,2) = use_p_index(use_sig_neuron,15);

pre_choice_correct_error(:,1) = use_neuron_index(use_sig_neuron,16); %Correct choice index
pre_choice_correct_error(:,2) = use_neuron_index(use_sig_neuron,17); %Error choice index
p_pre_choice_correct_error(:,1) = use_p_index(use_sig_neuron,16);
p_pre_choice_correct_error(:,2) = use_p_index(use_sig_neuron,17);

SameDif_preTone_correct_error(:,1) = use_neuron_index(use_sig_neuron,18);
SameDif_preTone_correct_error(:,2) = use_neuron_index(use_sig_neuron,19);
p_SameDif_preTone_correct_error(:,1) = use_p_index(use_sig_neuron,18);
p_SameDif_preTone_correct_error(:,2) = use_p_index(use_sig_neuron,19);


preCorrectTone_index_correct_correct_error(:,1) = use_neuron_index(use_sig_neuron,20);
preCorrectTone_index_correct_correct_error(:,2) = use_neuron_index(use_sig_neuron,21);
p_preCorrectTone_index_correct_correct_error(:,1) = use_p_index(use_sig_neuron,20);
p_preCorrectTone_index_correct_correct_error(:,2) = use_p_index(use_sig_neuron,21);

preCorrectChoice_index_correct_correct_error(:,1) = use_neuron_index(use_sig_neuron,22);
preCorrectChoice_index_correct_correct_error(:,2) = use_neuron_index(use_sig_neuron,23);
p_preCorrectChoice_index_correct_correct_error(:,1) = use_p_index(use_sig_neuron,22);
p_preCorrectChoice_index_correct_correct_error(:,2) = use_p_index(use_sig_neuron,23);


pre_same_pLcL_pLcH_index(:,1) = use_neuron_index(use_sig_neuron,24);
pre_same_pHcL_pHcH_index(:,1) = use_neuron_index(use_sig_neuron,25);
cur_same_pLcL_pHcL_index(:,1) = use_neuron_index(use_sig_neuron,26);
cur_same_pLcH_pHcH_index(:,1) = use_neuron_index(use_sig_neuron,27);

p_pre_same_pLcL_pLcH_index(:,1) = use_p_index(use_sig_neuron,24);
p_pre_same_pHcL_pHcH_index(:,1) = use_p_index(use_sig_neuron,25);
p_cur_same_pLcL_pHcL_index(:,1) = use_p_index(use_sig_neuron,26);
p_cur_same_pLcH_pHcH_index(:,1) = use_p_index(use_sig_neuron,27);


norm_spike = use_norm_spike(use_sig_neuron,:);
p_surprise = use_p_surprise(use_sig_neuron,:);

temp = norm_spike(:,2) - norm_spike(:,1); %right-left
neuron_choice_category = zeros(length(temp),1);
temp = temp > 0;
neuron_choice_category(temp) = 1; %right=1, left=0

%%%
%Remove NAN
test = [pre_sound, norm_spike, choice_correct_error,tone_correct_error,...
    pre_tone_correct_error,pre_choice_correct_error,SameDif_preTone_correct_error, p_SameDif_preTone_correct_error,...
    preCorrectTone_index_correct_correct_error,p_preCorrectTone_index_correct_correct_error,...
    preCorrectChoice_index_correct_correct_error,p_preCorrectChoice_index_correct_correct_error,...
    pre_same_pLcL_pLcH_index, p_pre_same_pLcL_pLcH_index,...
    pre_same_pHcL_pHcH_index,p_pre_same_pHcL_pHcH_index,...
    cur_same_pLcL_pHcL_index, p_cur_same_pLcL_pHcL_index,...
    cur_same_pLcH_pHcH_index,p_cur_same_pLcH_pHcH_index ];
test = mean(test,2);
test = isnan(test);
test_nan = find(test == 1);
pre_sound(test_nan,:) = [];
p_pre_sound(test_nan,:) = [];
norm_spike(test_nan,:) = [];
p_surprise(test_nan,:) = [];
choice_correct_error(test_nan,:) = [];
p_choice_correct_error(test_nan,:) = [];
neuron_choice_category(test_nan,:) = [];

tone_correct_error(test_nan,:) = []; %Correct choice index
p_tone_correct_error(test_nan,:) = [];

pre_tone_correct_error(test_nan,:) = []; %Correct choice index
p_pre_tone_correct_error(test_nan,:) = [];

pre_choice_correct_error(test_nan,:) = [];
p_pre_choice_correct_error(test_nan,:) = [];

SameDif_preTone_correct_error(test_nan,:) = [];
p_SameDif_preTone_correct_error(test_nan,:) = [];

preCorrectTone_index_correct_correct_error(test_nan,:) = [];
p_preCorrectTone_index_correct_correct_error(test_nan,:) = [];

preCorrectChoice_index_correct_correct_error(test_nan,:) = [];
p_preCorrectChoice_index_correct_correct_error(test_nan,:) = [];

pre_same_pLcL_pLcH_index(test_nan,:) = [];
p_pre_same_pLcL_pLcH_index(test_nan,:) = [];
pre_same_pHcL_pHcH_index(test_nan,:) = [];
p_pre_same_pHcL_pHcH_index(test_nan,:) = [];
cur_same_pLcL_pHcL_index(test_nan,:) = [];
p_cur_same_pLcL_pHcL_index(test_nan,:) = [];
cur_same_pLcH_pHcH_index(test_nan,:) = [];
p_cur_same_pLcH_pHcH_index(test_nan,:) = [];


return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sig_neuron = get_sig_neuron_time_window(p_task,new_p_thre,window)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
use_p_sound = min(p_task(:,window),[],2);
use_p_sound = -log10(use_p_sound);
sig_neuron = find(use_p_sound > new_p_thre);

return



