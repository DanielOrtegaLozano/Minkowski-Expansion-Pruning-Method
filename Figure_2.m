%% This code illustrates the MEP method by D. Ortega-Lozano & V. Makarov
%
% It is for the paper:
%
% Cognitive Maps and Path Planning in Minkowski Spacetime: 
% An Image-Dilation Approach
%
% This code simulates the drone's flight (Figure 4).  
%

clear

%% The method parameters
alpha = 3;       % steps ratio
h = 0.01;        % space step
tau = alpha * h; % time step
Tmax = 4;        % maximum look-ahead time
Kmax = floor(Tmax / tau);

%% Definition of obstacles
v = .7; % Maximum velocity of obstacles
x0 = 1; % initial position of the mobile obstacle

Mobs2 = inv(diag([.2, .1, .4]));
Mobs4 = inv(diag([.15, .2, .4]));
obstacle1 = @(t, G) vecnorm(G - [.1; 1; .6], Inf) <= .6;
obstacle2 = @(t, G) vecnorm(Mobs2 * (G - [.8; .6; .3]), Inf) <= 1;
obstacle3 = @(t, G) vecnorm(G - [.4; .2; .8], Inf) <= .2;
obstacle4 = @(t, G) vecnorm(Mobs4 * (G - [x0-v*t; .2; .3]), Inf) <= 1;
obstacle_union = @(t, G) obstacle1(t, G) | obstacle2(t, G) | obstacle3(t, G) | obstacle4(t, G);

[X, Y, Z] = meshgrid(0:h:1);
G = [X(:)'; Y(:)'; Z(:)'];
obstacle = @(t) reshape(obstacle_union(t, G), size(X));

% balls for dilation
S_alpha = makedisk(alpha);
beta = (1 + v / 2) * alpha + sqrt(3);
S_beta = makedisk(beta);

%% Step 1: Compute and thicken discretized obstacle in Minkowsky Space, Lambda
Lambda = false([size(X), Kmax]);
for k = 1:Kmax
    Lambda(:, :, :, k) = obstacle(tau * (k - 1/2));
end
Pure_Obstacle = Lambda;
Lambda = imdilate(Lambda, S_beta); % Dilation by beta along spatial dimensions only

%% Step 2: Evaluate the reachable spacetime domain, R
R = false([size(X), 1 + Kmax]);
R(1, 1, 1, 1) = true;
tic
for k = 1:Kmax
    R(:, :, :, k + 1) = imdilate(R(:, :, :, k), S_alpha) & ~Lambda(:, :, :, k);
end
toc

%% Step 3: Compute the Cognitive Map of the situation, T 
[ArrivalSteps, idx] = max(R, [], 4);
T = squeeze(idx - 1);
T(squeeze(~ArrivalSteps)) = Inf;

%% Find trajectory
% Target 
x = round(.9 /h); y = round(.98 /h); z = round(.10 /h);

Ki = T(x,y,z) + 1; % best time

B = false(size(X));
path = zeros(4, Ki);
path(:, Ki) = [x; y; z; Ki];
for k = Ki-1:-1:1
    Bk = B;
    Bk(x, y,  z) = true;
    Bk = imdilate(Bk, S_alpha);
    C = Bk & R(:, :, :, k);
    i = find(C, 1);
    [x, y, z] = ind2sub(size(X), i);
    path(:, k) = [x; y; z; k];
end

%% Visualization: Thickened obstacles
figure('Color','w') 
make_object_video(Pure_Obstacle)

%% Visualization: Building reachable set
figure('Color','w')
make_object_video(R)

%% Visualization: Trajectory in the configuration space
figure('Color','w') 
make_route_video(Pure_Obstacle, path)


%%%%%%% Auxilary functions %%%%%%%%%%%%%%%%





%% Make disk function
function S = makedisk(r)
    rceil = ceil(r);    
    [I, J, K] = meshgrid(-rceil:rceil, -rceil:rceil, -rceil:rceil);
    D = sqrt(I.^2 + J.^2 + K.^2);    
    S = double(D <= r);
end

%% Visualization tools
function view_spacetime_frame(A)
    A([1 end],:,:) = 0;
    A(:,[1 end],:) = 0;
    A(:,:,[1 end]) = 0;
        
    pch = patch(isosurface(double(A(:,:,:)), .5));
    set(pch, 'FaceColor', [0.5 1 0.5], ...
             'EdgeColor', 'none', ...
             'FaceAlpha', 0.8);   % transparency (0 = invisible, 1 = opaque)
    view(23.1000,   38.4000)
    camlight
    lighting phong
    xlim([1 size(A,1)])
    ylim([1 size(A,2)])
    zlim([1 size(A,3)])
    grid on
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function make_object_video(A)
    frames = round(linspace(1,size(A, 4),100));
    for i = frames
        clf
        view_spacetime_frame(A(:, :, :, i))
        title(['k = ' num2str(i)])
        drawnow
        pause(.1)
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function make_route_video(A, path)
    frames = round(linspace(1,size(A, 4),100));
    for i = frames
        path_past = path(4, :) <= i;
        truncated_path = path(:, path_past);
        clf
        view_spacetime_frame(A(:, :, :, i))
        hold on
        plot3(truncated_path(2,:), truncated_path(1,:), truncated_path(3,:), 'r-', LineWidth=4)
        plot3(truncated_path(2,end), truncated_path(1,end), truncated_path(3,end), 'r.', MarkerSize=25)
        plot3(path(2, end), path(1, end), path(3, end), 'b*', MarkerSize=10)
        title(['k = ' num2str(i)])
        drawnow
        pause(.1)
    end
end