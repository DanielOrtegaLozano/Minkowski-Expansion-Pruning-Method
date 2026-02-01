%% This code illustrates the MEP method by D. Ortega-Lozano & V. Makarov
%
% It is for the paper:
%
% Cognitive Maps and Path Planning in Minkowski Spacetime: 
% An Image-Dilation Approach
%
% This minimal code reproduces Figure 2.  
%

clear 

%% Definition of obstacles in Figure 2 (moving circles).
g1 = @(t, x) vecnorm(x - [.5; .5] ) < .3;
g2 = @(t, x) vecnorm(x - [1-.7*t; .25]) < .3;
g3 = @(t, x) vecnorm(x - [.4*t; .4+.3*t]) < .3;

g = @(t, x) g1(t, x) | g2(t, x) | g3(t, x);
v = .7; % Maximum velocity of obstacles

%% The algorithm's parameters
alpha = 10;  % ratio of time-space steps
h = 0.002;   % space discretization step
Tmax = 3;    % maximum look-ahead time

tau = alpha * h; % time step
Kmax = floor(Tmax / tau); % number of time steps
[X, Y] = meshgrid(0:h:1);
obstacle = @(t) reshape(g(t, [X(:)'; Y(:)']), size(X));

%% Step 1: Compute and thicken obstacle in the M-space
Lambda = false([size(X), Kmax]);
for k = 1:Kmax
    Lambda(:, :, k) = obstacle(tau * (k - 1/2));
end
beta = (1 + v / 2) * alpha + sqrt(2);
S_beta = makedisk(beta);
Lambda = imdilate(Lambda, S_beta);

%% Step 2: Evaluate the reachable set R
R = false([size(X), 1 + Kmax]);
R(1, 1, 1) = true;
S_alpha = makedisk(alpha);
for k = 1:Kmax
    R(:, :, k + 1) = imdilate(R(:, :, k), S_alpha) & ~Lambda(:, :, k);
end

%% Step 3: Compute the Cognitive Map
[ArrivalSteps,idx] = max(R,[],3);
T = idx-1; T(~ArrivalSteps) = Inf;

%% Find the fastest trajectory to a target and plot
target = [.95; .8];  % target's location

x = round(target(1) / h); y = round(target(2) / h); % end-point
Ki = T(x, y) + 1;    % best time
path = zeros(Ki,2);
path(Ki,:) = [x y];
for k = Ki-1:-1:1
    E = false(size(X));
    E(x, y) = true;
    E = imdilate(E, S_alpha) & R(:, :, k);
    [candidates_x, candidates_y] = ind2sub(size(X), find(E));
    [~, closest_candidate_ind] = min(vecnorm([candidates_x candidates_y]' - [x; y]));
    x = candidates_x(closest_candidate_ind); 
    y = candidates_y(closest_candidate_ind);
    path(k,:) = [x y];
end
path = path*h;

% Figure
figure('Color','w')
imagesc(0:h:1,0:h:1,T)
set(gca,'YDir','normal')
hold on
plot(path(:,1),path(:,2),'w','LineWidth',3)
plot(target(1),target(2),'o','Color','b',...
    'MarkerSize',10,'MarkerFaceColor','c')
colormap("turbo")
text(target(1)-0.05,target(2)+0.04,'target','FontSize',12,'Color','w')
axis square
title('The shortest trajectory in the cognitive map')

%% Make disk function
function S = makedisk(r)
    rceil = ceil(r);    
    [I, J] = meshgrid(-rceil:rceil, -rceil:rceil);
    D = sqrt(I.^2 + J.^2);    
    S = double(D <= r);
end
