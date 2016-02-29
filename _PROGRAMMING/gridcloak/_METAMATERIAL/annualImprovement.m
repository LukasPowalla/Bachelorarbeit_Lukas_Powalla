function [Q,ignoredPercentage] = annualImprovement(alphas,betas,P_values,varargin)
% Integrates yearly intensity transmitted to solar cell and finds the
% ratio Q for direct and/or diffuse irradiation.
% Q is the ratio of the integrated intensity with cloak vs. the integrated intensity
% without cloak. Q is calculated according to the formula
%
%   Q = integral(P(t)*T_top_substrate(t)*I_tot(t) dt) / integral(T_top_substrate(t)*I_tot(t) dt)
%
% The second output argument, ignoredPercentage, gives the percentage of
% the total annually deposited energy density deposited that comes in at
% incidence angles exceeding the range covered by alphas, betas. E.g. a
% value of 0.1 means 10% of the annual energy is not taken into account for
% Q computation, because there are no P_values for the respective angle
% combinations, i.e. the relative error of the result is 10%. To reduce
% this error, give P_values covering a wider angular range.
%
%  Arguments:
%   [Q,ignoredPercentage] = annualImprovement(alphas,betas,P_values): builds a
%       scatteredInterpolant from the vectors alphas,betas and P_values (same
%       length required) and computes Q.

addpath sunSimulator;
isacityClass = @(input) isa(input,'cityClass');
isaorientedIrradianceComputation = @(input) strcmpi(input,'global') || strcmpi(input,'dni');

p = inputParser;

p.addParameter('n1',1,@isnumeric);
p.addParameter('n2',4,@isnumeric);
p.addParameter('city',cityClass('UNSPECIFIED'),isacityClass);
p.addParameter('dayOfYearSweep',1:365,@isinteger);
p.addParameter('localTimeSweep',0:0.25:23.75,@isnumeric);
p.addParameter('cellAzimuth',6,@isnumeric);
p.addParameter('cellInclination',36,@isnumeric);
p.addParameter('cellRotation',90,@isnumeric);
p.addParameter('plots',false,@islogical);
p.addParameter('plotEveryNth',90,@isinteger);
p.addParameter('includeDiffusive',true,@islogical);
p.addParameter('includeDirect',true,@islogical);
p.addParameter('orientedIrradianceComputation','global',isaorientedIrradianceComputation);
p.addParameter('firstFigure',101,@isinteger);

p.parse(varargin{:});

n1 = p.Results.n1;
n2 = p.Results.n2;
city = p.Results.city;
dayOfYearSweep = p.Results.dayOfYearSweep;
localTimeSweep = p.Results.localTimeSweep;
cellAzimuth = p.Results.cellAzimuth;    % optimum for Karlsruhe: 6. 0 means cell facing south
cellInclination = p.Results.cellInclination; % optimum for Karlsruhe: 36.
cellRotation = p.Results.cellRotation;   % rotation of cell on roof. 0 means solar cell wires from bottom to top of roof
plots = p.Results.plots;
plotEveryNth = p.Results.plotEveryNth;
includeDiffusive = p.Results.includeDiffusive;  % include diffusive irradiation
includeDirect = p.Results.includeDirect;        % include direct irradiation
if strcmpi(p.Results.orientedIrradianceComputation,'dni') % compute irradiance on oriented plane using DNI or global irradiance
   fromDNI = true;
else
   fromDNI = false;
end

if strcmp(city.description,'UNSPECIFIED')
    city = cityClass('Karlsruhe',1,49.007,8.404,['sunSimulator' filesep 'solar_data' filesep 'karlsruhe'],fromDNI);
end

if includeDiffusive && fromDNI
       error('diffusive irradiation cannot be included when irradiance on oriented plane is computed using DNI');
end
firstFigure = p.Results.firstFigure;

P = scatteredInterpolant(alphas(:),betas(:),P_values(:),'linear','none');

function [f,a] = setupFigs(num)
   for k = 1:num
       f(k) = figure(firstFigure-1+k);
       clf(f(k));
       a(k) = gca;
       hold(a(k),'on');
   end
end

%% create figures
if plots
    [f,a] = setupFigs(8);
end

%% initialize variables
FF = 0.1;
solarCell = solarCellClass(cellAzimuth, cellInclination, cellRotation);
totSize = [length(localTimeSweep), length(dayOfYearSweep)];
mat = zeros(totSize);
sunHeight=mat; sunAzimuth=mat; beta=mat; alpha=mat; theta=mat;
totalIntensity=mat; totalIntensityDiff=mat; globalIntensity=mat;
cI = mat; cId = mat; bI = mat; bId = mat;

dt = mean(diff(localTimeSweep));

if includeDiffusive
    [TdiffCloak,TdiffBare,~] = FresnelDif(P,[],n1,n2,false);
else
    TdiffCloak = 0;
    TdiffBare = 0;
end

if includeDiffusive && ~includeDirect
    P = griddedInterpolant({[0 1],[0 1]},[0 0; 0 0],'linear','nearest');
end

n = 1;

%% transmittance from air to silicon
function T = T_top_substrate(angle)
    T = FresnelT(angle,[n1 n2],'none')';
end

%% loop through year

for dayOfYear = dayOfYearSweep
    [sunHeight(:,n), sunAzimuth(:,n)] = city.sunAngles(dayOfYear,localTimeSweep);
    [alpha(:,n), beta(:,n)] = solarCell.incidentAngles(city,dayOfYear,localTimeSweep);
    [totInt,totIntDiff,theta(:,n),~] = solarCell.totalIntensity(city,dayOfYear,localTimeSweep);
    if ~includeDiffusive
        totIntDiff(:,:) = 0;
    end
    if includeDiffusive && ~includeDirect
        totInt(:,:) = 0;
    end
    
    %% calculate intensities that are transmitted into the substrate
     % (factor (1-FF) cancels out at division of integrals (see labbook),
     % so 
    % bare intensity with direct irradiance (i.e. no cloak, no reference structure, with wire)
    bI(:,n) = (1-FF)*totInt'.*T_top_substrate(abs(theta(:,n)));
    
    % bare intensity with diffuse irradiance (i.e. no cloak, no reference structure, with wire)
    bId(:,n) = (1-FF)*totIntDiff'*TdiffBare;
    
    % intensity with direct irradiance and cloak in place
    cI(:,n) = (1-FF)*totInt'.*P(abs(alpha(:,n)),abs(beta(:,n))).*T_top_substrate(abs(theta(:,n))); 
    
    % intensity with diffuse irradiance and cloak in place
    cId(:,n) = (1-FF)*totIntDiff'*TdiffCloak;
    
    % intensities impinging onto the cell ("before" Fresnel transmittance,
    % without wire)
    totalIntensity(:,n) = totInt;
    totalIntensityDiff(:,n) = totIntDiff;
    globalIntensity(:,n) = totInt + totIntDiff;  
        
    n = n+1;
end

    function B = zeroNans(A)
        B = A;
        B(isnan(A)) = 0;
    end

%% calculate energy densities
globalEnergyDensity = globalIntensity*dt;
totalEnergyDensity = totalIntensity*dt;
totalEnergyDensityDiff = totalIntensityDiff*dt;
cloakEnergyDensity = (zeroNans(cI) + zeroNans(cId))*dt;
bareEnergyDensity = (zeroNans(bI) + zeroNans(bId))*dt;
dailyCloakEnergyDensity = sum(zeroNans(cloakEnergyDensity));   % in Wh/m^2/d
dailyBareEnergyDensity = sum(zeroNans(bareEnergyDensity));       % in Wh/m^2/d
accuCloakEnergyDensity = sum(zeroNans(dailyCloakEnergyDensity)); % in Wh/m^2/a
accuBareEnergyDensity = sum(zeroNans(dailyBareEnergyDensity));   % in Wh/m^2/a
Q = accuCloakEnergyDensity / accuBareEnergyDensity;

accuNannedEnergyDensity = sum(sum(zeroNans(totalEnergyDensity(isnan(cI))))) + sum(sum(zeroNans(totalEnergyDensityDiff(isnan(cId)))));
accuGlobalEnergyDensity = sum(sum(globalEnergyDensity));
ignoredPercentage = accuNannedEnergyDensity/accuGlobalEnergyDensity;


%% find annually deposited energy vs. (binned) incident angle
m = 1;
n = 1;
phis = 0:5:90;
binnedEnergyDensityBeta = zeros(length(phis),1);
binnedEnergyDensityAlpha = zeros(length(phis),1);

for d = dayOfYearSweep
    for t = localTimeSweep
        beta_now = beta(m,n);
        alpha_now = alpha(m,n);
        
        int_now = globalIntensity(m,n);
        
        [~,~,binBeta] = histcounts(beta_now,phis);
        [~,~,binAlpha] = histcounts(alpha_now,phis);
        
        if binBeta > 0 && binBeta <= length(phis)
            binnedEnergyDensityBeta(binBeta) = binnedEnergyDensityBeta(binBeta) + int_now*dt;
        end
        
        if binAlpha > 0 && binAlpha <= length(phis)
            binnedEnergyDensityAlpha(binAlpha) = binnedEnergyDensityAlpha(binAlpha) + int_now*dt;
        end
        
        m = m+1;
    end
    m = 1;
    n = n+1;
end

%% plotting
if plots
    dphi = mean(diff(phis));
    beta_lin = reshape(beta,1,numel(beta));
    alpha_lin = reshape(alpha,1,numel(alpha));
    globalIntensity_lin = reshape(globalIntensity,1,[]);
    alpha_hours = 365*24*histcounts(alpha,'BinWidth',5,'Normalization','probability');
    
    plot(a(1),sunAzimuth(:,1:plotEveryNth:end),sunHeight(:,1:plotEveryNth:end),'-o');
    plot(a(2),beta(:,1:plotEveryNth:end),alpha(:,1:plotEveryNth:end),'-o');
    plot(a(3),localTimeSweep,globalIntensity(:,1:plotEveryNth:end),'-o');
    plot(a(4),localTimeSweep,totalIntensityDiff(:,1:plotEveryNth:end),'-o');
    plot(a(5),dayOfYearSweep,dailyCloakEnergyDensity, dayOfYearSweep,dailyBareEnergyDensity);
    bar(a(6),[min(min(beta)):5:max(max(beta))],alpha_hours);
    scatter(a(7),beta_lin,globalIntensity_lin);
    scatter(a(7),alpha_lin,globalIntensity_lin);
    bar(a(8),phis-dphi*.5,1e-3*[binnedEnergyDensityBeta, binnedEnergyDensityAlpha]);
end


c = num2cell(dayOfYearSweep);
leg = {};
m = 1;
for j = 1:plotEveryNth:length(c)
    leg{m} = datestr(c{j}+366,'dd. mm.');
    m = m+1;
end

if plots
    legend(a(1),leg);
    xlabel(a(1),'Sun Azimuth (deg)');
    ylabel(a(1),'Sun Height (deg)');
    xlim(a(1),[0 360]);
    ylim(a(1),[0 90]);

    legend(a(2),leg);
    xlabel(a(2),'beta (deg)');
    ylabel(a(2),'alpha (deg)');

    legend(a(3),leg);
    xlabel(a(3),'local time (h)');
    ylabel(a(3),'global intensity on cell (no wire) (W/m^2)');

    legend(a(4),leg);
    xlabel(a(4),'local time (h)');
    ylabel(a(4),'diffusive intensity on cell (no wire) (W/m^2)');

    legend(a(5),{'with cloak','without cloak'});
    xlabel(a(5),'day of year');
    ylabel(a(5),'energy density on cell with wire (Wh/m^2/d)'); 
    
    xlabel(a(6),'alpha (deg)');
    ylabel(a(6),'duration spent in alpha during full year (h/a)');
   
    legend(a(7),{'beta','alpha'});
    xlabel(a(7),'angle (deg)');
    ylabel(a(7),'intensity on blank cell (W/m^2)');
    
    legend(a(8),{'beta','alpha'});
    xlabel(a(8),'angle (deg)');
    ylabel(a(8),'energy density deposited at angle during full year (kWh/m^2/a)');
end

end