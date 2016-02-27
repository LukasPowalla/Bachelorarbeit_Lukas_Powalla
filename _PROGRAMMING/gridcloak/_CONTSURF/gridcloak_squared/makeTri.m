function [P,N,tr] = makeTri(X0,Y0,Z0)
%% triangulate meshgrid
% INPUTS:
%  X0,Y0: x- and y-coordinates of 2D meshgrid
%  Z0: heights associated to X0 and Y0
%
% OUTPUTS:
%  P: center points of triangles
%  N: surface normals of triangles
%  tr: triangulation object

doPlot = false;

s = size(X0);
P = [reshape(X0,[],1) reshape(Y0,[],1) reshape(Z0,[],1)];

tri = [];

for n = 0:s(1)-2
    for m = 1:s(1):s(1)*(s(2)-1)
        mat = [m m+s(1) m+1;
            m+s(1) m+s(1)+1 m+1] + n*ones(2,3);
        tri = [tri; mat];
    end
end

tr = triangulation(tri,P);
P = tr.incenter([1:size(tr,1)]');
N = tr.faceNormal;

if doPlot
    figure(111)
    trisurf(tr);
    hold on
    scatter3(P(:,1),P(:,2),P(:,3),'r*')
    line([P(:,1)'; P(:,1)'+N(:,1)'], [P(:,2)'; P(:,2)'+N(:,2)'], [P(:,3)'; P(:,3)'+N(:,3)']);
    axis equal
end