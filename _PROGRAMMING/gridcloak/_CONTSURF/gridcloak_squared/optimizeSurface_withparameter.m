function [vopt] = optimizeSurface_withparameter(varargin)
ng=4;
%% input parsing
p = inputParser;

p.addParameter('a',1,@isnumeric);
p.addParameter('b',1,@isnumeric);
p.addParameter('y0x',5,@isnumeric);
p.addParameter('y0y',5,@isnumeric);                 
p.addParameter('rx1',5,@isnumeric);        
p.addParameter('rx2',50,@isnumeric);            
p.addParameter('ry1',5,@isnumeric);  
p.addParameter('ry2',50,@isnumeric);


p.parse(varargin{:});

v(1)=p.Results.a;
v(2)=p.Results.b;
v(3)=p.Results.y0x;
v(4)=p.Results.y0y;


%% fitness and plotting functions
    function [fitness,v] = getFitness(v)
        a =v(1);
        b =v(2);
        y0x =v(3);
        y0y =v(4);
        fitness=-getQofFreeformwithParameter('a',a,'b',b,'y0x',y0x,'y0y',y0y);
    end

    function state = plotSurf(options,state,flag)
        try
           vcur = state.Population(state.Score==min(state.Score),:);
           bar(vcur(1,:));
           set(gca,'XTickLabel',{'a','b','y0x','y0y'})
           xlabel('');
           ylabel('value');
           title(sprintf('improvement = %2.3f',-min(state.Score)));
        end
    end

%% set optimization options
options = gaoptimset(...
    'PlotFcns',{@plotSurf},... % @gaplotbestindiv
    'FitnessLimit',-1.5,...
    'Generations',inf,...
    'PlotInterval',1,...
    'PopulationSize',50,...
    'CreationFcn',@gacreationuniform ...
    );
%% optimize!
v_LB=[ 0 0 0 0];
v_UB=[ 5 5 10 10];

[vopt,fitness] = ga(@getFitness,ng,[],[],[],[],v_LB,v_UB,[],options);
display(fitness)

toc
end