%%   Parameter

amin=0;
amax=7;
bmin=0;
bmax=6;
y0xmin=30;
y0xmax=80;
y0ymin=30;
y0ymax=80;
an=10;
bn=10;
Yn=2;
%%
a=linspace(amin,amax,an);
b=linspace(bmin,bmax,bn);
y0x_s=linspace(y0xmin,y0xmax,Yn);
y0y_s=linspace(y0ymin,y0ymax,Yn);
improvement=zeros(an,bn);%,Yn,Yn
for i=1:an
    for k=1:bn
         improvement(i,k)=getQofFreeformwithParameter('a',a(i),'b',b(k));
%         for l=1:Yn
%             for m=1:Yn
%                    improvement(i,k,l,m)=getQofFreeformwithParameter('a',a(i),'b',b(k),'y0x',y0x_s(l),'y0y',y0y_s(m));
%             end
%         end
    end
end
maximum=max(improvement(:))
[A,B]=meshgrid(linspace(amin,amax,an),linspace(bmin,bmax,bn));
%% plotting
surf(A,B,improvement);
title(sprintf('max. improvement = %2.3f',max(improvement(:))));
xlabel('a')
ylabel('b')
zlabel('improvement');
shading interp;
colorbar;