%% make main and supplementary figures
% This is the summary program for making the main and supplementary figures
% Please change the root dir here and functions for data dir 

root_dir = '/Volumes/Extreme SSD';

%% Figure 1 & Figure S1

Figure1HMM_psycho_full_model_Gauss_20240710
Figure1E_S1C_HMM_all_session_analysis_20230314_plot
Figure1F_G_H_nBackTrial_learningphase_20240420_aki(root_dir)

close all
% 
%% Figure 2 & Figure S1
% 
Figure2D_top_ExampleSessionRL_20240115(root_dir)
Figure2D_bottom_ExampleSessionState_20240115(root_dir)
Figure2E_S1B_HMM_all_session_20240111_model_plot_each_session
Figure2F_G_H_I_S1F_HMM_all_session_20240420_model_plot(root_dir)

close all

%% Figure 3 & Figure S1

Figure3A_auto2_HMM_behave_analysis_block3_plot
Figure3B_auto2_HMM_behave_analysis_block3_plot
Figure3C_D_nBackTrial_20240420_overtrain
Figure3E_F_G_S1I_HMM_all_session_analysis_20240423_overtrain

close all

%% Figure 4

Figure4B_S3B_Rand_trial_activity_20240417
Figure4D_HMM_20240522_regress_CV_compare_shuffle_repeat
Figure4D_HMM_20240522_regress_CV_compare_shuffle_alternate
Figure4E_F_GLM_CV_proportions(root_dir)

close all

%% Figure 5 & Figure S5A

Figure5B_Example_PriorTrace_20230919(root_dir)
Figure5C_S5A_process_HMM_20240711_before_sound_depth
Figure5D_E_process_HMM_20240703_surprise_before_sound_depth
Figure5G_process_HMM_20240710_compare_RL

close all

%% Figure 6 & Figure S5B,C

Figure6A_S5B_20240710_scatter_before_sound_depth
Figure6B_S5C_20240710_trace_before_sound_depth
Figure6C_6D_20240710_all_regrion_before_sound_depth
Figure6E_ME_ephys_regress_20240628(root_dir)

close all

%% Figure 7 & Figure S6A,B,C

Figure7B_ExampleNeuron_DuringSound(root_dir)
Figure7C_S6A_HMM_20231129_index_at_sound_depth
Figure7D_S6B_HMM_20231129_scatter_at_sound_depth
Figure7E_S6C_HMM_20230510_trace_at_sound_depth_prefer
Figure7F_G_HMM_20231129_all_region_at_sound_depth

close all

%% Figure 8 & Figure S7

Figure8B_ChoiceTrace_20230919(root_dir)
Figure8C_S7A_HMM_20230504_index_choice_ver2_depth_correct
Figure8D_S7B_HMM_20230504_scatter_ver2_depth_correct
Figure8E_S7C_HMM_20230504_trace_choice_depth_prefer
Figure8F_G_HMM_20230504_all_region_depth_correct
Figure8H_ME_ephys_regress_step2_20240628_correct

close all

%% Figure S3

FigureS3D_Neg_neuron_ave_trace2_depth_step2
FigureS3E_GLM_CV_proportions_decreasing(root_dir)
FigureS4F_HMM_20240521_regress_CV_delta_proportion

close all

%% Figure S4

FigureS4B_top_MotionEnergy_example_session_20240521(root_dir)
FigureS4B_bottom_MotionEnergy_CorrectError_step2_20240521
FigureS4C_MotionenergyRegression_20240521
FigureS4D_MotionEnergyCV_summary

close all

%% Figure S5D,E

FigureS5D_OFC_scatter_current_outcome_before_sound
FigureS5E_CurrentOutcome_20240703_surprise_before_sound_depth

close all

%% Figure S6D,E

FigureS6D_OFC_scatter_current_outcome_during_sound
FigureS6E_CurrentOutcome_20240705_surprise_at_sound_all

close all

%% Figure S8

FigureS8A_B_LeftPanel_HMM_20230504_scatter_choice_error
FigureS8A_B_RightPanel_HMM_20230504_trace_choice_error
FigureS8C_D_HMM_20230504_all_regions_depth_error
FigureS8E_ME_regress_step2_20240628_error

close all
