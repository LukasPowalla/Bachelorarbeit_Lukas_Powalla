addpath('export_figure');

set(gcf,'Color','w');
set(gca,'Color','w');

pos = get(gcf,'Position');

clear xlim ylim
title('')
axis equal
xlim([-80 80]);
ylim([-80 80]);


screenRes = 177;

width = 10;
height = 10;
width=width/2.54/1.844;
height=height/2.54/1.844;

set(gcf,'Position',[pos(1:2) width*screenRes height*screenRes]);
set(gca,'FontSize',10);

export_fig sd.png -png -r600 %-nocrop