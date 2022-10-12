function gillespie()
%% Rate constants
p.kR = 0.016;      
p.kP = 0.016;
p.gR = 0.0101;                        
p.gP = 0.02536;
%% Initial state
tspan = [0, 10000]; 
x0    = [0, 0];     
%% Specify reaction network
pfun = @propensities_2state;
stoich_matrix = [ 1  0    %MAKING Y2-Y3
                 -1  0    %y2-y3 with mir 10b
                  0  1    %decay-Y3/initiator
                  0 -1 ]; %using initiator in hcr
%% Run simulation
[t,x] = directMethod(stoich_matrix, pfun, tspan, x0, p);
%% Plot time course
figure();
stairs(t,x); set(gca,'XLim',tspan);
xlabel('time (s)');
ylabel('molecules');
legend({'Y2-Y3','Initiator'});
end
function a = propensities_2state(x, p)
% Return reaction propensities given current state x
mRNA    = x(1);
protein = x(2);
a = [p.kR;            %transcription
     p.kP*mRNA;       %translation
     p.gR*mRNA;       %mRNA decay
     p.gP*protein;]   %protein decay
end
function [ t, x ] = directMethod( stoich_matrix, propensity_fcn, tspan, x0 ,params)
%% Initialize
num_species = size(stoich_matrix,2);
T = zeros(1000000, 1);
X = zeros(1000000, num_species);
T(1)     = tspan(1);
X(1,:)   = x0;
rxn_count = 1;
%% MAIN LOOP
while T(rxn_count) < tspan(2)        
    % Calculate reaction propensities
    a = propensity_fcn(X(rxn_count,:), params);
    
    % Sample earliest time-to-fire (tau)
    a0 = sum(a);
    r = rand(1,2);
    tau = (1/a0)*log(1/r(1));
    
    % Sample identity of earliest reaction channel to fire (mu)
    mu=1; s=a(1); r0=r(2)*a0;
    while s < r0
       mu = mu + 1;
       s = s + a(mu);
    end
    if rxn_count + 1 > 1000000
        t = T(1:rxn_count);
        x = X(1:rxn_count,:);
        return;
    end
    
    % Update time and carry out reaction mu
    T(rxn_count+1)   = T(rxn_count)   + tau;
    X(rxn_count+1,:) = X(rxn_count,:) + stoich_matrix(mu,:);    
    rxn_count = rxn_count + 1;
    
end  
% Return simulation time course
t = T(1:rxn_count);
x = X(1:rxn_count,:);
if t(end) > tspan(2)
    t(end) = tspan(2);
    x(end,:) = X(rxn_count-1,:);
end    
end