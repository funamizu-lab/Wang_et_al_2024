
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
function [analysis_dir,depth_def] = repeat_PPC_20240427

depth_def = 1;

analysis_dir = {
    '/Volumes/Extreme SSD/Zigzag_repeat_Ephys/E_phys/repeat/W03/PPC/2021-09-23_11-31-16_unknown_errors_W03_PPC_L_2nd_ReAnalyze'
    '/Volumes/Extreme SSD/Zigzag_repeat_Ephys/E_phys/repeat/W03/PPC/2021-09-24_11-31-29_done_W03_PPC_L_3rd'
    '/Volumes/Extreme SSD/Zigzag_repeat_Ephys/E_phys/repeat/W03/PPC/2021-10-13_10-02-32_done_W03_PPC_R_1st_LICK_IS_NOT_GOOD'
    '/Volumes/Extreme SSD/Zigzag_repeat_Ephys/E_phys/repeat/W03/PPC/2021-10-14_09-58-17_done_W03_PPC_R_2nd'
    '/Volumes/Extreme SSD/Zigzag_repeat_Ephys/E_phys/repeat/W03/PPC/2021-10-15_10-01-00_unknown_errors_W03_PPC_R_3rd_ReAnalyze'

    '/Volumes/Extreme SSD/Zigzag_repeat_Ephys/E_phys/repeat/W06/PPC/2021-11-17_09-38-29_W06_PPC_L_1st_errorMakingSpikeFiles'
    '/Volumes/Extreme SSD/Zigzag_repeat_Ephys/E_phys/repeat/W06/PPC/2021-11-18_09-26-10_W06_PPC_L_2nd'
    '/Volumes/Extreme SSD/Zigzag_repeat_Ephys/E_phys/repeat/W06/PPC/2021-11-19_09-35-31_W06_PPC_L_3rd'
    '/Volumes/Extreme SSD/Zigzag_repeat_Ephys/E_phys/repeat/W06/PPC/2021-12-09_10-20-55_W06_PPC_R_2nd'
    '/Volumes/Extreme SSD/Zigzag_repeat_Ephys/E_phys/repeat/W06/PPC/2021-12-10_10-10-18_W06_PPC_R_3rd'
    
    '/Volumes/Extreme SSD/Zigzag_repeat_Ephys/E_phys/repeat/W22/PPC/2022-03-16_10-50-06_W22_PPC_L_1st'
    '/Volumes/Extreme SSD/Zigzag_repeat_Ephys/E_phys/repeat/W22/PPC/2022-03-17_13-17-12_W22_PPC_L_2nd'
    '/Volumes/Extreme SSD/Zigzag_repeat_Ephys/E_phys/repeat/W22/PPC/2022-03-18_09-22-08_W22_PPC_L_3rd'

    '/Volumes/Extreme SSD/Zigzag_repeat_Ephys/E_phys/repeat/W31/PPC/2022-11-16_14-10-42_W31_PPC_L_1st_AfterLED'
    '/Volumes/Extreme SSD/Zigzag_repeat_Ephys/E_phys/repeat/W31/PPC/2022-11-17_15-02-34_W31_PPC_L_2nd_AfterLED'
    '/Volumes/Extreme SSD/Zigzag_repeat_Ephys/E_phys/repeat/W31/PPC/2022-11-18_14-05-01_W31_PPC_L_3rd_AfterLED'
    '/Volumes/Extreme SSD/Zigzag_repeat_Ephys/E_phys/repeat/W31/PPC/2022-12-07_16-03-52_W31_PPC_R_1st_AfterLED'
    '/Volumes/Extreme SSD/Zigzag_repeat_Ephys/E_phys/repeat/W31/PPC/2022-12-08_15-41-46_W31_PPC_R_2nd_AfterLED'
    };

