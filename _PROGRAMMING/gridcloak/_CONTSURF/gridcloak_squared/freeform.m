function Z = freeform(R1,R2,y0,X)
%returns the numerical 1D solurion with the parameter of the solution R1,
%R2, y0 refering to the x-values X
n=1.5;

Xend = R2;

dxMin = 0.001;
dxMax = 0.5;
dxRamp = 10; % number of steps from dxMin to dxMax

do_plot = false;

max_height = [];


    dx = linspace(dxMin,dxMax,dxRamp);
    s = cumsum(dx);
    if s(end)>R2
        dx = dx(1,find(s<=R2,1,'last'));
    else
        dx = [dx repmat(dxMax,size(s(end):dxMax:R2))];
    end
    x0 = [0 cumsum(dx(1:end-1))];

    
    %% calculate approximate z(x) (small angles approximation)
    z0small = sqrt(y0^2+R1/(1-1/n)*(2*x0-x0.^2/R2));
    g0small = atand(diff(z0small)./diff(x0));
    x0small = x0+0.5*dx;
    x0small = x0small(1:end-1);
    
    %% calculate g and z for variable height (exact numeric solution)
    imp_fun = @(g,I,g0,x0,i) R1*(1-abs(x0(i+1))/R2)/tand(g-asind(sind(g)/n))-I(i)-tand(g0(i))*(x0(i+1)-x0(i));
        
    g0 = zeros(1,length(x0));
    I = zeros(1,length(x0));
    
    g0_first_fun = @(g) R1./(tand(g-asind(sind(g)/n)))-y0;
    g0(1) = fzero(g0_first_fun,g0small(1));
       
    if g0(1)<90    
        I(1) = y0;

        ifun = @(g) imp_fun(g,I,g0,x0,1);
        g0(2) = fzero(ifun,g0small(2));

        for i = 2:length(x0)-2

            I(i) = I(i-1) + tand(g0(i-1))*(x0(i)-x0(i-1));
            ifun = @(g) imp_fun(g,I,g0,x0,i);
            g0(i+1) = fzero(ifun,g0small(i+1));

            I(i+1) = I(i) + tand(g0(i+1))*(x0(i+1)-x0(i));  
        end

        I(i+2) = I(i+1);
        g0(i+2) = 0;



        %% plot g(x)
        if do_plot
            figure;

            plot(x0/R2,g0,x0small/R2,g0small,'LineWidth',1);%,x0,g0approx);
            xlabel('\it{x}/\it{R_2}');
            ylabel('\beta (deg)');
            legend({'full numerical solution','small angles approximation'});
        end


        %% plot z(x)
        if do_plot
            figure;

            plot([x0 Xend]/R2,[I I(end)]/y0,[x0 Xend]/R2,[z0small z0small(end)]/y0,'LineWidth',1);%,[x0 Xend],[z0approx' z0approx(end)]);
            xlabel('\it{x}/\it{R_2}');
            ylabel('\it{y(x)}/\it{y(0)}');
            legend({'numerical solution','analytical approximation'},'Location','SouthEast');
        end
    else
%         warning('angle above 90°'); % DEBUG
        disp(y0);
    end
    
    Z = interp1([x0 Xend],[I I(end)],X(:));
    
    Z = flip(Z);
    
    Z = reshape(Z,size(X));
    
    
end