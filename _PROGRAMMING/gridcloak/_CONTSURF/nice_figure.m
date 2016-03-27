addpath('export_figure');

set(gcf,'Color','w');
set(gca,'Color','w');

pos = get(gcf,'Position');

clear xlim ylim

% xlim([-50 50]);
% ylim([-50 50]);
% axis equal

screenRes = 177;

width = 20;
height = 10;
width=width/2.54/1.844;
height=height/2.54/1.844;

set(gcf,'Position',[pos(1:2) width*screenRes height*screenRes]);
set(gca,'FontSize',10);

export_fig test.png -png -r600 -nocrop