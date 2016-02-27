function [ Schnittpunkte,Normale ] = calculateSPofRayswithSurface( tr,tr_Points1,tr_Points2,tr_Points3,orig,v_orig )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

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
sp_list=NaN(m_max,k_max,3);


%number of triangle is m
%number of ray is k

for k=1:k_max
    for m=1:m_max
        t=(fn(m,1)*(tr_Points1(m,1)-orig(k,1))+fn(m,2)*(tr_Points1(m,2)-orig(k,2))+fn(m,3)*(tr_Points1(m,3)-orig(k,3)))/((fn(m,1)*v_orig(k,1))+(fn(m,2)*v_orig(k,2))+(fn(m,3)*v_orig(k,3)));
        sp_list(m,k,1)=(orig(k,1)+t*v_orig(k,1));
        sp_list(m,k,2)=(orig(k,2)+t*v_orig(k,2));
        sp_list(m,k,3)=(orig(k,3)+t*v_orig(k,3));
        a=tr_Points2(m,1)-tr_Points1(m,1);
        b=tr_Points3(m,1)-tr_Points1(m,1);
        c=sp_list(m,k,1)-tr_Points1(m,1);
        d=tr_Points2(m,2)-tr_Points1(m,2);
        e=tr_Points3(m,2)-tr_Points1(m,2);
        f=sp_list(m,k,2)-tr_Points1(m,2);
        g=tr_Points2(m,3)-tr_Points1(m,3);
        h=tr_Points3(m,3)-tr_Points1(m,3);
        i=sp_list(m,k,3)-tr_Points1(m,3);
        t1=(e*c-b*f)/(a*e-d*b);
        t2=(a*f-d*c)/(a*e-d*b);
        if isnan(t1)||isnan(t2)
            if (a==0)||(b==0)
                t1=(h*f-e*i)/(d*h-e*g);
                t2=(i*d-g*f)/(d*h-e*g);
            elseif (d==0)||(e==0)
                t1=(c*h-b*i)/(a*h-g*b);
                t2=(a*i-g*c)/(a*h-g*b);
            end
            
        end
        if t>(10^(-10))
            if (t1>=0&&t2>=0&&t1+t2<=1)
                t_list(m,k)=t;
            end
        end
    end
    m=find((min(t_list(:,k))==t_list(:,k)),1);
     if (~isempty(m))&&(t_list(m,k)>0)
        Schnittpunkte(k,:)=sp_list(m,k,:);
        Normale(k,:)=fn(m,:);
     end

end

end
