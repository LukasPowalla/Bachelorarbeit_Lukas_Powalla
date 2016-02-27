
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

amin=2;
amax=2.5;
h0min=20;
h0max=100;
Y0min=20;
Y0max=90;
an=5;
h0n=5;
Y0n=5;

[A,H0,Y0] = meshgrid(linspace(amin,amax,an),linspace(h0min,h0max,h0n),linspace(Y0min,Y0max,Y0n));

Av=A(:);
H0v=H0(:);
Y0v=Y0(:);


fit = [];

for i=1:length(Av)   
 fit(i) = getFitness_freeform(Av(i),H0v(i),Y0v(i)); 
end
   
fit = reshape(fit,size(A));
minimum=zeros([Y0n,1]);
aktuellesminimum=min(min(fit(:,:,1)));

for k=1:Y0n
    minimum(k)=min(min(fit(:,:,k)));
    if aktuellesminimum>=min(min(fit(:,:,k)));
        aktuellesminimum=min(min(fit(:,:,k)));
        y0i=k;
    end
end
surf(A(:,:,y0i),H0(:,:,y0i),fit(:,:,y0i));
title(sprintf('y0 = %2.1f\nmin. fitness = %2.1f',Y0(1,1,y0i),min(min(fit(:,:,y0i)))));
xlabel('a')
ylabel('h0')
zlabel('fitness');

