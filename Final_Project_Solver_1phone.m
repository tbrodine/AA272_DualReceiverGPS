% Matthew Chmiel and Taylor Brodine
% AA 272 Final Project
clc, clear, close all

% Read in the log data, preprocessed by ProcessGnssMeasScript
log = readtable('groundtruth_l1.csv');

%% divide data into chunks based on RxTime_s
start_ids = 1;
for i = 2:length(log.X)
    if log.Rx_time(i) == log.Rx_time(i-1)
        continue
    else
        start_ids = [start_ids,i];
    end
end
fprintf('%0.3f seconds logged\n',log.Rx_time(end)-log.Rx_time(1));
fprintf('%d time steps logged\n',length(start_ids));

%% loop through all time instances
% initialize solution vectors
pos = zeros(4,length(start_ids));
pos_nl = zeros(4,length(start_ids));
for i = 1:length(start_ids)
    if i~=length(start_ids)
        ids = start_ids(i):start_ids(i+1)-1;
    else
        ids = start_ids(i):size(log,1);
    end
    
    x0 = [0;0;0;0];
    % Solve using N-R
%     pos(:,i) = solve_pos(x0,log.X(ids),log.Y(ids),log.Z(ids),...
%         log.B(ids),log.rho(ids)); % rho is called Pseudo_m in gnss_log
    % Solve using N-R with Weighted LS
    pos(:,i) = solve_pos_WLS(x0,log.X(ids),log.Y(ids),log.Z(ids),...
        log.B(ids),log.rho(ids),log.sigma_rho(ids));
    
    % Solve using nonlinear cost function:
%     fun = @(x) cost(x,log.X(ids),log.Y(ids),log.Z(ids),...
%         log.B(ids),log.rho(ids));
%     options = optimoptions(@fmincon,'display','off');
%     pos_nl(:,i) = fmincon(fun,x0,[],[],[],[],[],[],[],options);    
end

truth = [mean(pos(1,:));mean(pos(2,:));mean(pos(3,:))];
lla_pos = ecef2lla(pos(1:3,:)')';
lla_truth = [mean(lla_pos(1,:));mean(lla_pos(2,:));mean(lla_pos(3,:))];

figure
plot3(pos(1,:),pos(2,:),pos(3,:),'-b')
hold on
plot3(truth(1),truth(2),truth(3),'xr','MarkerSize',20,'linewidth',4)
xlabel('x position [m]')
ylabel('y position [m]')
zlabel('z position [m]')
title('ECEF Trajectory of Positions')
grid on
legend({'WLS','Mean'})
axis equal

figure
plot(lla_pos(2,:),lla_pos(1,:),'-b')
hold on
plot(lla_truth(2),lla_truth(1),'xr','MarkerSize',20,'linewidth',4)
xlabel('Longitude [deg W]')
ylabel('Latitude [deg N]')
title('LLA Trajectory of Positions')
grid on
legend({'WLS','Mean'})
axis equal

save('truth_l1.mat','truth')
fprintf('Mean Values from Phone 0 Truth run:\n')
fprintf('X = %0.6e m\n',truth(1))
fprintf('Y = %0.6e m\n',truth(2))
fprintf('Z = %0.6e m\n\n',truth(3))
fprintf('Latitude = %0.6e deg N\n',lla_truth(1))
fprintf('Longitude = %0.6e deg W\n',lla_truth(2))
fprintf('Altitude = %0.6e m\n\n',lla_truth(3))

%% Supporting Functions:

% cost function for NL optimization solver:
function [c] = cost(xa,Xa,Ya,Za,Ba,rho_a)%,xb,Xb,Yb,Zb,Bb,rho_b)
% ca = [];
% for i = 1:length(Xa)
%     ca = [ca;rho_a(i)-sqrt((xa(1)-Xa(i))^2+(xa(2)-Ya(i))^2+(xa(3)-Za(i))^2)...
%         -xa(4)+Ba(i)];
% end
ca = norm(rho_a-get_expected_pseudoranges(xa,Xa,Ya,Za,Ba));
c = norm(ca);
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
