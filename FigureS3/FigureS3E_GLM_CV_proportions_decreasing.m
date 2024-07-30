function FigureS3E_GLM_CV_proportions_decreasing(root_dir)

RepBefore = [0	0	0	0	0.0323	0.0645	0.3548
0	0.0667	0	0.1333	0	0	0.1333
0.0357	0.0357	0.0357	0.0357	0	0.1429	0.3214
0	0.5	0	0	0	0	0.5];

RepDuring = [0.0556	0.0833	0.0278	0.0278	0	0.0833	0.2778
0	0.0588	0	0	0	0.0588	0.2353
0.0504	0.2773	0.0168	0.0042	0.0084	0.1092	0.2227
0	0.25	0	0	0	0.25	0 ];
RepChoice = [0.1127	0.0563	0.3239	0.0282	0.0141	0.0563	0.1549
0.0727	0.0364	0.4364	0.0545	0	0.1273	0.2364
0.1478	0.1773	0.4138	0.0197	0.0345	0.0887	0.2611
0.1515	0.0606	0.5758	0	0	0.2121	0.303 ];

ZigBefore = [0.1329	0.0153	0.0204	0.1244	0.0443	0.1261	0.4702
0.1458	0.0042	0.0208	0.0792	0.0208	0.075	0.4667
0.0657	0.0058	0.0117	0.0993	0.0409	0.2073	0.4175
0.1435	0.0093	0.0046	0.1157	0.0417	0.1111	0.5648
0.0472	0.0063	0.0031	0.0629	0.0189	0.1164	0.3899
0.0553	0	0	0.0632	0.0158	0.1395	0.3605 ];

ZigDuring = [0.1837	0.0424	0.0251	0.1068	0.0188	0.0879	0.3344
0.2215	0.0163	0.013	0.0847	0.0228	0.0717	0.4169
0.0947	0.387	0.0245	0.0462	0.016	0.0674	0.2591
0.2305	0.0532	0.0142	0.1064	0.0426	0.078	0.4929
0.0963	0.0248	0.0155	0.0404	0.0217	0.1056	0.3292
0.1121	0.0117	0.0117	0.0444	0.0117	0.1098	0.2991 ];

ZigChoice = [0.3991	0.2042	0.4714	0.0945	0.014	0.0782	0.2625
0.2449	0.1224	0.398	0.0969	0.0153	0.0306	0.1531
0.2531	0.1566	0.4547	0.0457	0.0145	0.0341	0.1958
0.3562	0.1644	0.5114	0.0776	0.0228	0.0822	0.3744
0.3421	0.2018	0.4123	0.0497	0.0175	0.0526	0.3129
0.2989	0.125	0.5326	0.0543	0.0054	0.0489	0.2065 ];


[green, purple] = getColor(root_dir);


figure
heatmap(RepBefore,'ColorLimits',[0 0.75],'Colormap',green)
title('Before Sound (Repeating condition)')
figure
heatmap(RepDuring,'ColorLimits',[0 0.75],'Colormap',green)
title('During Sound (Repeating condition)')
figure
heatmap(RepChoice,'ColorLimits',[0 0.75],'Colormap',green)
title('During Choice (Repeating condition)')

figure
heatmap(ZigBefore,'ColorLimits',[0 0.75],'Colormap',purple)
title('Before Sound (Alternating condition)')
figure
heatmap(ZigDuring,'ColorLimits',[0 0.75],'Colormap',purple)
title('During Sound (Alternating condition)')
figure
heatmap(ZigChoice,'ColorLimits',[0 0.75],'Colormap',purple)
title('During Choice (Alternating condition)')
return


function [green, purple] = getColor(root_dir)

path1 = strcat(root_dir, '/IntermediateFunctions/green_color_map.txt');
path2 = strcat(root_dir, '/IntermediateFunctions/purple_color_map.txt');
data1=importdata(path1);
data2=importdata(path2);
green=data1(:,1:3);
purple=data2(:,1:3);

return










