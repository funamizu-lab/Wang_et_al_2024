

function Figure3B_auto2_HMM_behave_analysis_block3_plot

data_rep = load('repeat_sabun_psycho_all_each_session');
plot_sabun_each_session(data_rep)

data_zig = load('zigzag_sabun_psycho_all_each_session');
plot_sabun_each_session(data_zig)



return


function plot_sabun_each_session(data)

opt_low1 = data.opt_low1;
opt_low2 = data.opt_low2;
opt_low3 = data.opt_low3;
opt_high1 = data.opt_high1;
opt_high2 = data.opt_high2;
opt_high3 = data.opt_high3;


evi_x = 0:0.01:1;
figure

plot_mean_se_moto_x_axis(opt_low1, evi_x, [51 153 51]./255,0) %Green
hold on
plot_mean_se_moto_x_axis(opt_high1, evi_x, [255 191 17]./255,0) %Yellow
hold on
plot_mean_se_moto_x_axis(opt_low2, evi_x, [0 204 255]./255,0) %Cyan
hold on
plot_mean_se_moto_x_axis(opt_high2, evi_x, [232 108 76]./255,0) %dark red
hold on
plot_mean_se_moto_x_axis(opt_low3, evi_x, [0 0 1],0) %Blue
hold on
plot_mean_se_moto_x_axis(opt_high3, evi_x, [1 0 0],0) %red
set(gca,'xlim',[-0.1 1.1], 'ylim', [0 1])

return



