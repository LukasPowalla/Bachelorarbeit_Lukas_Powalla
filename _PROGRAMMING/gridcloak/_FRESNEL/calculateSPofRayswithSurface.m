function [ Schnittpunkte,Normale ] = calculateSPofRayswithSurface( tr,tr_Points1,tr_Points2,tr_Points3,orig,v_orig )
%Intersection of the rays defined through orig and v_orig with tr
%  we calculate the 'closest' intersection from the rays with
%  the triangulatio object.

fn=faceNormal(tr);
m_max=size(fn,1);
k_max=size(orig,1);
for i=1:m_max
    if (fn(i,3)<0==1)
        fn(i,:)=-fn(i,:);
    end
end
%% intersection from one ray with one triangular

Normale=NaN(k_max,3);
Schnittpunkte=NaN(k_max,3);
t_list=NaN(m_max,k_max);

%number of triangle is m
%number of ray is k
% t describes the parameter for the Intersection from the ray with the
% plane of the m-th triangle. In order to know, whether the rays hits the
% m-th triangle, we have to calculate the parameters t1 and t2. This can be
% seen as follows:
% SP=t1*P_12+t2*P_13            where:          SP= P_start+t*v_start
% t has to be positive and t1>=0&&t2>=0&&t1+t2<=1 in order to have a
% intersection. Then, pick the smallest t-value because the first
% inersection with respect to t is the one, which is hit first in reality.
% In the end, we return the "Schnittpunkt" and the "Normale" referring to
% the Intersection point. 

% a=nan(m_max,1);
% d=nan(m_max,1);
% g=nan(m_max,1);

% b=nan(m_max,1);
% c=nan(m_max,1);
% e=nan(m_max,1);
% f=nan(m_max,1);
% h=nan(m_max,1);
% i=nan(m_max,1);
% u=nan(m_max,1);
% v=nan(m_max,1);
% w=nan(m_max,1);
    b=tr_Points2(:,1)-tr_Points1(:,1);
    c=tr_Points3(:,1)-tr_Points1(:,1);
    e=tr_Points2(:,2)-tr_Points1(:,2);
    f=tr_Points3(:,2)-tr_Points1(:,2);
    h=tr_Points2(:,3)-tr_Points1(:,3);
    i=tr_Points3(:,3)-tr_Points1(:,3);
    
for k=1:k_max
    a=-v_orig(k,1);
    d=-v_orig(k,2);
    g=-v_orig(k,3);
    u=orig(k,1)-tr_Points1(:,1);
    v=orig(k,2)-tr_Points1(:,2);
    w=orig(k,3)-tr_Points1(:,3);
    
    det=a.*e.*i-a.*f.*h-b.*d.*i+b.*f.*g+c.*d.*h-c.*e.*g;
    t=(b.*f.*w-b.*i.*v-c.*e.*w+c.*h.*v+e.*i.*u-f.*h.*u)./det;
    t1=(-a.*f.*w+a.*i.*v+c.*d.*w-c.*g.*v-d.*i.*u+f.*g.*u)./det;
    t2=(a.*e.*w-a.*h.*v-b.*d.*w+b.*g.*v+d.*h.*u-e.*g.*u)./det;
    t0=find(t>1e-10 & t1>=0 & t2>=0 & t1+t2<=1);
    t_list(t0,k)=t(t0);
end

M = nan(1,k_max);
T=min(t_list(:,:));
A=repmat(T,size(t_list,1),1);
B=(A)==t_list(:,:);
[row,col]=find(B); %,1,'first'
M(col) = row;
Schnittpunkte(:,:)=(orig(:,:)+v_orig(:,:).*repmat(T,3,1)');
Normale(~isnan(M),:)=fn(M(~isnan(M)),:);
end
