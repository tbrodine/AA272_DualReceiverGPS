% Matthew Chmiel and Taylor Brodine
% AA 272 Final Project
clc, clear, close all

% Read in the log data, preprocessed by ProcessGnssMeasScript
% Uncomment one of the blocks below, depending on what data to analyze

% % % TEST 1
% log0 = readtable('phone0_1m_l1.csv');
% log1_raw = readtable('phone1_1m_l1.csv');
% global set_distance
% set_distance = 1; % m
% moving = 0;
% load truth_l1.mat truth
% titl = 'Test 1 - ';

% % % TEST 2
% log0 = readtable('phone0_5m_l1.csv');
% log1_raw = readtable('phone1_5m_l1.csv');
% global set_distance
% set_distance = 5; % m
% moving = 0;
% load truth_l1.mat truth
% titl = 'Test 2 - ';

% % % TEST 3
% log0 = readtable('phone0_1m_l2.csv');
% log1_raw = readtable('phone1_1m_l2.csv');
% global set_distance
% set_distance = 1; % m
% moving = 0;
% load truth_l2.mat truth
% titl = 'Test 3 - ';

% % % TEST 4
% log0 = readtable('phone0_5m_l2.csv');
% log1_raw = readtable('phone1_5m_l2.csv');
% global set_distance
% set_distance = 5; % m
% moving = 0;
% load truth_l2.mat truth
% titl = 'Test 4 - ';

% % % TEST 5
% log0 = readtable('right_35_4ft.csv');
% log1_raw = readtable('left_35_4ft.csv');
% global set_distance
% set_distance = 4*0.3048; % m
% moving = 1;
% truth = 35*0.44704; % m/s
% titl = 'Test 5 - ';

% % % TEST 6
% log1_raw = readtable('right_35_4ft.csv');
% log0 = readtable('left_35_4ft.csv');
% global set_distance
% set_distance = 4*0.3048; % m
% moving = 1;
% truth = 35*0.44704; % m/s
% titl = 'Test 6 - ';

% keep
if moving==0
lla_truth = ecef2lla(truth')';
end

%% divide data into chunks based on RxTime_s
start_ids0 = 1;
for i = 2:length(log0.X)
    if log0.Rx_time(i) == log0.Rx_time(i-1)
        continue
    else
        start_ids0 = [start_ids0,i];
    end
end
% fprintf('Phone 0:\n%0.3f seconds logged\n',log0.Rx_time(end)-log0.Rx_time(1));
% fprintf('%d time steps logged\n\n',length(start_ids0));

%% interpolate phone1 data onto phone0 time instances - for moving
if moving==1
    % loop through each SVID used in phone1
    svids1 = unique(log1_raw.SVID);
    t0 = unique(log0.Rx_time);
    t1 = unique(log1_raw.Rx_time);
    log1 = log1_raw;
    for i = 1:length(svids1)
        idx = log1_raw.SVID==svids1(i);
        if sum(idx)~=length(t0)
            t0i = zeros(sum(idx),1);
            t1i = log1.Rx_time(idx);
            for j = 1:sum(idx)
                [~,idx_t0] = min(abs(t0-t1i(j)));
                t0i(j) = t0(idx_t0);
            end
            t0i = sort(t0i);
        else
            t0i=t0;
        end
        log1.X(idx) = interp1(log1_raw.Rx_time(idx),log1_raw.X(idx),t0i,'linear','extrap');
        log1.Y(idx) = interp1(log1_raw.Rx_time(idx),log1_raw.Y(idx),t0i,'linear','extrap');
        log1.Z(idx) = interp1(log1_raw.Rx_time(idx),log1_raw.Z(idx),t0i,'linear','extrap');
        log1.B(idx) = interp1(log1_raw.Rx_time(idx),log1_raw.B(idx),t0i,'linear','extrap');
        log1.rho(idx) = interp1(log1_raw.Rx_time(idx),log1_raw.rho(idx),t0i,'linear','extrap');
        log1.Rx_time(idx) = t0i;
    end
    [~,order] = sort(log1.Rx_time);
    log1 = log1(order,:);
else
    log1 = log1_raw;
end

%% divide phone1 data into chunks based on time

start_ids1 = 1;
for i = 2:length(log1.X)
    if log1.Rx_time(i) == log1.Rx_time(i-1)
        continue
    else
        start_ids1 = [start_ids1,i];
    end
end
% fprintf('Phone 1:\n%0.3f seconds logged\n',log1.Rx_time(end)-log1.Rx_time(1));
% fprintf('%d time steps logged\n\n',length(start_ids1));

%% loop through all time instances for phone0
% initialize solution vectors
pos0_alone = zeros(4,length(start_ids0));
% pos_nl0 = zeros(4,length(start_ids0));
for i = 1:length(start_ids0)
    if i~=length(start_ids0)
        ids = start_ids0(i):start_ids0(i+1)-1;
    else
        ids = start_ids0(i):size(log0,1);
    end
    
    x0 = [0;0;0;0];
    % Solve using N-R
%     pos0_alone(:,i) = solve_pos(x0,log0.X(ids),log0.Y(ids),log0.Z(ids),...
%         log0.B(ids),log0.rho(ids)); % rho is called Pseudo_m in gnss_log
    % Solve using N-R
    pos0_alone(:,i) = solve_pos_WLS(x0,log0.X(ids),log0.Y(ids),log0.Z(ids),...
        log0.B(ids),log0.rho(ids),log0.sigma_rho(ids)); % rho is called Pseudo_m in gnss_log
end
    pos0_alone_mean = [mean(pos0_alone(1,:));mean(pos0_alone(2,:));mean(pos0_alone(3,:))];
lla_pos0_alone = ecef2lla(pos0_alone(1:3,:)')';
lla_pos0_alone_mean = [mean(lla_pos0_alone(1,:));mean(lla_pos0_alone(2,:));mean(lla_pos0_alone(3,:))];

figure(1)
plot3(pos0_alone(1,:),pos0_alone(2,:),pos0_alone(3,:),'-b')
hold on
if moving==0
   plot3(pos0_alone_mean(1),pos0_alone_mean(2),pos0_alone_mean(3),'xb','MarkerSize',20,'linewidth',4)
end

figure(2)
plot(lla_pos0_alone(2,:),lla_pos0_alone(1,:),'-b')
hold on
if moving==0
   plot(lla_pos0_alone_mean(2),lla_pos0_alone_mean(1),'xb','MarkerSize',20,'linewidth',4)
end

%% phone 1 alone as the IC
pos1_alone = zeros(4,length(start_ids1));
% pos_nl0 = zeros(4,length(start_ids0));
for i = 1:length(start_ids1)
    if i~=length(start_ids1)
        ids = start_ids1(i):start_ids1(i+1)-1;
    else
        ids = start_ids1(i):size(log1,1);
    end
    
    x0 = [0;0;0;0];
    % Solve using N-R
%     pos0_alone(:,i) = solve_pos(x0,log0.X(ids),log0.Y(ids),log0.Z(ids),...
%         log0.B(ids),log0.rho(ids)); % rho is called Pseudo_m in gnss_log
    % Solve using N-R
    pos1_alone(:,i) = solve_pos_WLS(x0,log1.X(ids),log1.Y(ids),log1.Z(ids),...
        log1.B(ids),log1.rho(ids),log1.sigma_rho(ids)); % rho is called Pseudo_m in gnss_log
end

%% compute phone0 position but using NL opt. w/ Phone 2
time_steps = min([length(start_ids0),length(start_ids1)]);
% initialize solution vectors
pos0_NL = zeros(4,time_steps);
pos1_NL = zeros(4,time_steps);
for i = 1:time_steps
    if i~=length(start_ids0)
        ids0 = start_ids0(i):start_ids0(i+1)-1;
    else
        ids0 = start_ids0(i):size(log0,1);
    end
    if i~=length(start_ids1)
        ids1 = start_ids1(i):start_ids1(i+1)-1;
    else
        ids1 = start_ids1(i):size(log1,1);
    end
    
    x0 = [pos0_alone(:,i),pos1_alone(:,i)];
    
    % Solve using nonlinear cost function:
    fun = @(x) cost(x,...
        log0.X(ids0),log0.Y(ids0),log0.Z(ids0),...
        log0.B(ids0),log0.rho(ids0),...
        log1.X(ids1),log1.Y(ids1),log1.Z(ids1),...
        log1.B(ids1),log1.rho(ids1));
    options = optimoptions(@fmincon,'display','none',...
        'MaxFunctionEvaluations',10000,...
        'MaxIterations',5000,...
        'StepTolerance',1e-20,...
        'ConstraintTolerance',1e-3,...
        'OptimalityTolerance',1e-16,...
        'Algorithm','interior-point');
    nonlncon = @dist_constraint;
    [x_opt,fval] = fmincon(fun,x0,[],[],[],[],[],[],nonlncon,options);
%     fval
    pos0_NL(:,i) = x_opt(:,1);
    pos1_NL(:,i) = x_opt(:,2);
end
pos0_NL_mean = [mean(pos0_NL(1,:));mean(pos0_NL(2,:));mean(pos0_NL(3,:))];
lla_pos0_NL = ecef2lla(pos0_NL(1:3,:)')';
lla_pos0_NL_mean = [mean(lla_pos0_NL(1,:));mean(lla_pos0_NL(2,:));mean(lla_pos0_NL(3,:))];

figure(1)
plot3(pos0_NL(1,:),pos0_NL(2,:),pos0_NL(3,:),'-r')
if moving==0
   plot3(pos0_NL_mean(1),pos0_NL_mean(2),pos0_NL_mean(3),'xr','MarkerSize',20,'linewidth',4)
   plot3(truth(1),truth(2),truth(3),'xk','MarkerSize',20,'linewidth',4)
end
xlabel('x position [m]')
ylabel('y position [m]')
zlabel('z position [m]')
title([titl, 'ECEF Trajectory of Phone 0 Position'])
grid on
if(moving==0)
    legend({'WLS, Phone 0 Alone','WLS Mean','DUAL, Both Phones','DUAL Mean','Truth Mean'},'location','best')
else
legend({'WLS, Phone 0 Alone','DUAL, Both Phones'},'location','best')
end
axis equal

figure(2)
plot(lla_pos0_NL(2,:),lla_pos0_NL(1,:),'-r')
if moving==0
   plot(lla_pos0_NL_mean(2),lla_pos0_NL_mean(1),'xr','MarkerSize',20,'linewidth',4)
   plot(lla_truth(2),lla_truth(1),'xk','MarkerSize',20,'linewidth',4)
end
xlabel('Longitude [deg W]')
ylabel('Latitude [deg N]')
title([titl, 'LLA Trajectory of Phone 0 Position'])
grid on
if(moving==0)
    legend({'WLS, Phone 0 Alone','WLS Mean','DUAL, Both Phones','DUAL Mean','Truth Mean'},'location','best')
else
legend({'WLS, Phone 0 Alone','DUAL, Both Phones'},'location','best')
end
axis equal

if moving==0
%     fprintf('Mean Values from Phone 0 Alone:\n')
%     fprintf('X = %0.6e m\n',pos0_alone_mean(1))
%     fprintf('Y = %0.6e m\n',pos0_alone_mean(2))
%     fprintf('Z = %0.6e m\n\n',pos0_alone_mean(3))
%     fprintf('Mean Values from Phone 0 using Phone 1 data:\n')
%     fprintf('X = %0.6e m\n',pos0_NL_mean(1))
%     fprintf('Y = %0.6e m\n',pos0_NL_mean(2))
%     fprintf('Z = %0.6e m\n\n',pos0_NL_mean(3))
    fprintf('Phone 0 Alone XYZ distance from truth:\n')
    fprintf('Error = %0.3f m\n\n',norm(truth-pos0_alone_mean))
    fprintf('Phone 0 using Phone 1 XYZ distance from truth:\n')
    fprintf('Error = %0.3f m\n\n',norm(truth-pos0_NL_mean))
    fprintf('Phone 0 Alone Lat/Long distance from truth:\n')
    fprintf('Error = %0.3e deg\n\n',norm(lla_truth(1:2)-lla_pos0_alone_mean(1:2)))
    fprintf('Phone 0 using Phone 1 Lat/Long distance from truth:\n')
    fprintf('Error = %0.3e deg\n\n',norm(lla_truth(1:2)-lla_pos0_NL_mean(1:2)))
    
else
    vel0_alone = compute_vel_XYZ(pos0_alone);
    vel0_NL = compute_vel_XYZ(pos0_NL);
    lla_vel0_alone = compute_vel_LLA(lla_pos0_alone);
    lla_vel0_NL = compute_vel_LLA(lla_pos0_NL);
    fprintf('True Velocity from Car/Waze = %0.4f m/s\n\n',truth)
    fprintf('Phone 0 Alone XYZ Velocity:\n')
    fprintf('V = %0.4f m/s\n\n',mean(vel0_alone))
    fprintf('Phone 0 using Phone 1 XYZ Velocity:\n')
    fprintf('V = %0.4f m/s\n\n',mean(vel0_NL))
    fprintf('Phone 0 Alone Lat/Long Velocity:\n')
    fprintf('V = %0.4f m/s\n\n',mean(lla_vel0_alone))
    fprintf('Phone 0 using Phone 1 Lat/Long Velocity:\n')
    fprintf('V = %0.4f m/s\n\n',mean(lla_vel0_NL))
end

%% Supporting Functions:

% cost function for NL optimization solver:
function [c] = cost(x,Xa,Ya,Za,Ba,rho_a,Xb,Yb,Zb,Bb,rho_b)
xa = x(:,1);
xb = x(:,2);
ra = norm(rho_a-get_expected_pseudoranges(xa,Xa,Ya,Za,Ba));
rb = norm(rho_b-get_expected_pseudoranges(xb,Xb,Yb,Zb,Bb));
c = ra+rb;%sqrt(sum(ra.^2+rb.^2));
% ra = rho_a-get_expected_pseudoranges(xa,Xa,Ya,Za,Ba);
% rb = rho_b-get_expected_pseudoranges(xb,Xb,Yb,Zb,Bb);
% c = sqrt(sum(ra.^2)+sum(rb.^2));
end

% nonlinear constraint function:
function [c,ceq] = dist_constraint(x)
global set_distance
xa = x(:,1);
xb = x(:,2);
c = [];
ceq = norm(xa(1:3)-xb(1:3))-set_distance;
end

% Newton-Raphson method
function [x_vec] = solve_pos(x0,X,Y,Z,B,prange)
x_vec = x0;
del_y = [100;100;100;100];
while norm(del_y)>1e-3
    G = get_geomery_matrix(x_vec,X,Y,Z);
    prange_theo = get_expected_pseudoranges(x_vec,X,Y,Z,B);
    del_rho = prange-prange_theo;
    del_y = (((G')*G)\(G'))*del_rho;
    x_vec = x_vec+del_y;
end
end

% Construct the Geometry Matrix:
function [G] = get_geomery_matrix(x_est,X,Y,Z)
G = zeros(length(X),4);
for i = 1:length(X)
    G(i,:) = [-1*(([X(i);Y(i);Z(i)]-x_est(1:3))/...
        norm([X(i);Y(i);Z(i)]-x_est(1:3)))',1];
end
end

% Compute the expected p-ranges
function [rho] = get_expected_pseudoranges(x_est,X,Y,Z,B)
rho = zeros(length(X),1);
for i = 1:length(X)
    rho(i) = sqrt(([X(i);Y(i);Z(i)]-x_est(1:3))'*...
        ([X(i);Y(i);Z(i)]-x_est(1:3)))+x_est(4)-B(i);
end
end

function [vel] = compute_vel_XYZ(pos)
    vel = zeros(length(pos),1);
    for i = 1:length(pos)-1
        vel(i) = norm(pos(1:3,i)-pos(1:3,i+1));
    end
end

function [vel] = compute_vel_LLA(pos)
    r_E = 6371.000e3; % m
    c_E = 2*pi*r_E; % m
    vel = zeros(length(pos),1);
    for i = 1:length(pos)-1
        delta_lat = abs(pos(1,i)-pos(1,i+1))*pi/180;
        delta_long = abs(pos(2,i)-pos(2,i+1))*pi/180;
        dist_lat = c_E*delta_lat/2/pi; % m
        dist_long = (2*pi*r_E*cosd(pos(1,i)))*delta_long/2/pi; % m
        vel(i) = sqrt(dist_lat^2+dist_long^2);
    end
end

% Newton-Raphson method
function [x_vec] = solve_pos_WLS(x0,X,Y,Z,B,prange,sigma_rho)
x_vec = x0;
W = diag(1./sigma_rho);
del_y = [100;100;100;100];
while norm(del_y)>1e-3
    G = get_geomery_matrix(x_vec,X,Y,Z);
    prange_theo = get_expected_pseudoranges(x_vec,X,Y,Z,B);
    del_rho = prange-prange_theo;
    del_y = (((G')*W*G)\(G'))*W*del_rho;
    x_vec = x_vec+del_y;
end
end
