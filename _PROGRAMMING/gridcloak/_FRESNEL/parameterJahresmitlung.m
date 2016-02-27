
%%   Parameter
h0min=100;
h0max=350;
h0n=10;

%%
h0=linspace(h0min,h0max,h0n);

improvement=zeros(h0n);
for i=1:h0n
         improvement(i)=getQwithParameter('h0',h0(i));
end
maximum=max(improvement(:))
%% plotting
%[A,B]=meshgrid(linspace(amin,amax,an),linspace(bmin,bmax,bn));
% surf(A,B,improvement);
% title(sprintf('max. improvement = %2.3f',max(improvement(:))));
% xlabel('a')
% ylabel('b')
% zlabel('improvement');
% shading interp;
% colorbar;