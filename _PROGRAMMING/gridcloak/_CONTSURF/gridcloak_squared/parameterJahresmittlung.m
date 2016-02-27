%%   Parameter

amin=0;
amax=3;
bmin=0;
bmax=4;
y0xmin=3;
y0xmax=20;
y0ymin=3;
y0ymax=20;
an=10;
bn=10;
Yn=10;
%%
a=linspace(amin,amax,an);
b=linspace(bmin,bmax,bn);
y0x_s=linspace(y0xmin,y0xmax,Yn);
y0y_s=linspace(y0ymin,y0ymax,Yn);
improvement=zeros(an,bn);%,Yn,Yn
for i=3:an
    for k=1:6
%          improvement(i,k)=getQofFreeformwithParameter('a',a(i),'b',b(k));
        for l=1:10
            for m=1:10
                   improvement(i,k,l,m)=getQofFreeformwithParameter('a',a(i),'b',b(k),'y0x',y0x_s(l),'y0y',y0y_s(m)); 
            end
        end
        i
        k
    end
end
maximum=max(improvement(:))
% [A,B]=meshgrid(linspace(amin,amax,an),linspace(bmin,bmax,bn));
% %% plotting
% surf(A,B,improvement);
% title(sprintf('max. improvement = %2.3f',max(improvement(:))));
% xlabel('a')
% ylabel('b')
% zlabel('improvement');
% shading interp;
% colorbar;