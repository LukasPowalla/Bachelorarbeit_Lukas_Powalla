function [Zopt,X,Y] = optimizeSurface(varargin)
% optimizes surface with generic algorithm by degrees of freedom given by the number of Gridpoints
%% input parsing
p = inputParser;

p.addParameter('n1',1,@isnumeric);              % starting ref. index
p.addParameter('n2',1.5,@isnumeric);            % end ref. index
p.addParameter('h0',50,@isnumeric);             % minimum surface height
p.addParameter('finalDist',1,@isnumeric);       % stop optimization when average distance is smaller
p.addParameter('xlim',500,@isnumeric);          % x-size of unit cell
p.addParameter('ylim',50,@isnumeric);           % y-size of unit cell
p.addParameter('nx',50,@isnumeric);             % number of facets along x
p.addParameter('ny',5,@isnumeric);              % number of facets along y
p.addParameter('sx',0.9,@isnumeric);            % x-scaling factor
p.addParameter('sy',0.9,@isnumeric);            % y-scaling factor
p.addParameter('safetyFactor',0.9,@isnumeric);  % additional safety factor
p.addParameter('X0',[],@isnumeric);             % x-coordinates of initial surface
p.addParameter('Y0',[],@isnumeric);             % y-coordinates of initial surface
p.addParameter('Z0',[],@isnumeric);             % z-coordinates of initial surface

p.parse(varargin{:});

tic;
n1 = p.Results.n1;
n2 = p.Results.n2;
h0 = p.Results.h0;
finalDist = p.Results.finalDist;
xlim = p.Results.xlim;
ylim = p.Results.ylim;

nx = p.Results.nx;
ny = p.Results.ny;
sx = p.Results.sx;
sy = p.Results.sy;
safetyFactor = p.Results.safetyFactor;

X0 = p.Results.X0;
Y0 = p.Results.Y0;
Z0 = p.Results.Z0;

[X,Y] = meshgrid(linspace(0,xlim,nx),linspace(0,ylim,ny));


%% set starting surface
if isempty(X0)
%     Z0 = -h0*ones(size(X));
    ampl = 300;
    Z0 = ampl*(X.^2.+Y.^2)/max(X(:).^2+Y(:).^2)-h0-ampl;
else
    Z0 = interp2(X0/max(X0(:))*max(X(:)),Y0/max(Y0(:))*max(Y(:)),Z0,X,Y);
end

%% triangulate surface
P0 = makeTri(X,Y,Z0);

%% get sizes
np = size(P0,1);
ng = numel(Z0);

%% fitness and plotting functions
    function [fitness,Xhit,Yhit,tr,numOnContactGrid] = getFitness(Z)
        Z = reshape(Z,size(Z0));
        
        [P,N,tr] = makeTri(X,Y,Z);
        alpha=00;
        v0 = [sin(alpha) 0 cos(alpha)];
        v0 = repmat(v0,np,1);
        Xhit = nan*ones(np,1);
        Yhit = nan*ones(np,1);

        v1=Snell(v0,N,n1,n2);
        A=repmat([-1 0; 0 -1; 0 0],1,1,np);
        A(:,3,:)=reshape(v1',3,1,[]);

        Xr = P(:,1);
        Yr = P(:,2);
        Zr = P(:,3);
        B = [-Xr, -Yr, -Zr];
        B = reshape(B',3,1,[]);

        for m = 1:size(B,3)
            res = A(:,:,m)\B(:,:,m);
            Xhit(m) = res(1);
            Yhit(m) = res(2);
        end

        idealXr = safetyFactor*sx*Xr;
        idealYr = safetyFactor*sy*Yr;
        
        sqDist = (Xhit-idealXr).^2 + (Yhit-idealYr).^2;
        sqDistMod = sqDist .* (ones(size(sqDist)) + sqrt((Xhit.^2 + Yhit.^2)/(xlim^2+ylim^2)));

        fitness = sum(sum(sqDistMod));
        numOnContactGrid = length(find(Xhit>sx*xlim | Yhit>sy*ylim | Xhit<0 | Yhit<0));
        fitness = fitness * (numOnContactGrid+1);
    end

    function state = plotSurf(options,state,flag)
        try
           Zcur = state.Population(state.Score==min(state.Score),:);
           Zcur = reshape(Zcur,size(Z0));
           surf(X,Y,Zcur);
           colormap default;
        end
    end

    function state = plotDist(options,state,flag)
        try
           Zcur = state.Population(state.Score==min(state.Score),:);
           Zcur = Zcur(1,:);
           Zcur = reshape(Zcur,size(Z0));
           [~,Xhit,Yhit,tr,numOnContactGrid] = getFitness(Zcur);
           P = tr.incenter([1:size(tr,1)]');
           Xr = P(:,1);
           Yr = P(:,2);
           dist = reshape(sqrt((Xhit-safetyFactor*sx*Xr).^2+(Yhit-safetyFactor*sy*Yr).^2),[],1);
           patch('Faces',tr.ConnectivityList,'Vertices',tr.Points(:,1:2),'FaceVertexCData',dist);
           shading flat;
%            pcolor(X,Y,reshape(sqrt((Xhit-sx*reshape(X,ng,1)).^2+(Yhit-sy*reshape(Y,ng,1)).^2),size(X))); 
           title(sprintf('Distance from design: average %2.1f\n%d/%d rays hitting grid',mean(dist),numOnContactGrid,np));
           colorbar;
           colormap jet;
        end
    end

%% set optimization options
options = gaoptimset(...
    'PlotFcns',{@gaplotbestf, @plotSurf, @plotDist},...
    'FitnessLimit',finalDist.^2*np,...
    'Generations',inf,...
    'PlotInterval',1,...
    'PopulationSize',500,...
    'PopInitRange',[-2*h0;-h0],...
    'InitialPopulation',reshape(Z0,1,[]),...
    'CreationFcn',@gacreationuniform ...
    );

%% linear constraints
aeq = zeros(1,numel(Z0));
aeq(end) = 1;
beq = -h0;

%% optimize!
[Zopt,fitness] = ga(@getFitness,ng,[],[],[],[],-Inf,repmat(-h0,ng,1),[],options);

Zopt = reshape(Zopt,size(Z0));

toc
end