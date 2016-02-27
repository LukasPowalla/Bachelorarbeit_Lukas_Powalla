
% 
%%   Parameter

amin=1.8;
amax=2;
bmin=1.8;
bmax=2.1;
Y0min=70;
Y0max=80;
an=5;
bn=5;
Y0n=5;
%% spacing eqidistant Parameters
[A,B,Y0] = meshgrid(linspace(amin,amax,an),linspace(bmin,bmax,bn),linspace(Y0min,Y0max,Y0n));

Av=A(:);
Bv=B(:);
Y0v=Y0(:);

%% calculating fitness for all given Parameters 
%You can choose between the Euclidic mixing (sqrt(f(x)^2+g(x)^2)) of "analytical"/"numerical"-1D solution 
%and 

fit = [];

for i=1:length(Av)
    
 %fit(i) = getFitness(Av(i),Bv(i),Y0v(i));   
 fit(i) = getFitness_freeform(Av(i),Bv(i),Y0v(i));
 %fit(i) = getFitness_linear(Av(i),Bv(i),Y0v(i));
    
end
   
fit = reshape(fit,size(A));
%% searching best Y0
minimum=zeros([Y0n,1]);
aktuellesminimum=min(min(fit(:,:,1)));

for k=1:Y0n
    minimum(k)=min(min(fit(:,:,k)));
    if aktuellesminimum>=min(min(fit(:,:,k)));
        aktuellesminimum=min(min(fit(:,:,k)));
        y0i=k;
    end
end

%minI = find(fit==min(fit),1,'first');


%% plotting
surf(A(:,:,y0i),B(:,:,y0i),fit(:,:,y0i));
title(sprintf('y0 = %2.1f\nmin. fitness = %2.1f',Y0(1,1,y0i),min(min(fit(:,:,y0i)))));
xlabel('a')
ylabel('b')
zlabel('fitness');

