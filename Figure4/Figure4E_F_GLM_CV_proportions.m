function Figure4E_F_GLM_CV_proportions(root_dir)

RepBefore = [0.0143	0.0143	0.0286	0.1333	0.0714	0.2476	0.4571
0.25	0.25	0	0	0	0.75	0.25
0.0514	0.014	0.028	0.0607	0.0234	0.2243	0.5701
0.1304	0.0435	0.0435	0.1304	0.0435	0.2174	0.5652];

RepDuring = [0.2059	0.1429	0.0336	0.0882	0.0672	0.2227	0.4832
0	0	0	0	0	0.4286	0.2857
0.1087	0.4269	0.0632	0.0316	0.0198	0.1285	0.3063
0.2647	0.1176	0	0.0294	0.0294	0.1765	0.6176 ];

RepChoice = [0.3007	0.2108	0.5637	0.0572	0.0098	0.1503	0.2745
0.1598	0.0888	0.4556	0.0414	0.0059	0.1006	0.3432
0.222	0.1951	0.5852	0.0538	0.0123	0.139	0.278
0.2568	0.1577	0.4414	0.0315	0.009	0.1171	0.4505];

ZigBefore = [0.2236	0.022	0.0268	0.2016	0.0772	0.1827	0.6394
0.0339	0	0	0.0847	0.0508	0.0678	0.6949
0.0894	0.0073	0.0088	0.0968	0.0396	0.1378	0.5689
0.1141	0	0.0201	0.2282	0.0268	0.1409	0.7383
0.1046	0.0084	0.0209	0.1883	0.0377	0.2008	0.431
0.1176	0	0.0168	0.1261	0.0672	0.1681	0.4958 ];

ZigDuring = [0.4206	0.112	0.0575	0.1664	0.0484	0.1498	0.5809
0.2	0.0824	0.0118	0.0588	0.0118	0.0588	0.6118
0.1977	0.4687	0.0585	0.0544	0.0148	0.0651	0.3072
0.3147	0.066	0.0102	0.132	0.066	0.0711	0.7005
0.3221	0.1049	0.0262	0.1273	0.0337	0.1685	0.4831
0.2345	0.0759	0.0345	0.1103	0.0138	0.1655	0.4069 ];

ZigChoice = [0.4707	0.2667	0.6263	0.098	0.0131	0.1	0.3606
0.2513	0.098	0.5804	0.1156	0.0151	0.0352	0.2437
0.3551	0.2531	0.6586	0.0562	0.0078	0.0477	0.2217
0.4183	0.1662	0.6218	0.0831	0.0258	0.0401	0.3897
0.4971	0.2766	0.6132	0.0774	0.0155	0.1064	0.441
0.3736	0.1705	0.6372	0.0977	0.0155	0.1783	0.4031 ];




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










