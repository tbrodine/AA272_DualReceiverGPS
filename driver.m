

%% load dat
clc; clear;
prFileName = 'right_35_4ft.txt';
dirName = 'C:\MATLAB\Grad School\AA272\project';
[dat, QOIs] = PreProcessGNSS(dirName,prFileName);
writetable(QOIs,['data/' prFileName(1:end-3) 'csv'])

% dat = dat(dat.Cn0DbHz>15,:);

%% Solve

x0 = [0;0;0;0];
idx = 1;
tstep = 0;

% solve with their xyz and our rho
% while dat.RxTime_s(idx) <= dat.RxTime_s(end)-1 % set data pulling condition
%     tstep = tstep + 1;
%     % store current time step starting index
%     idx_start = idx; 
%     % check rows until time step changes
%     time_inst = dat.RxTime_s(idx_start);
%     while dat.RxTime_s(idx) == time_inst
% 
%         idx = idx+1;
%     end
%     % store current time step ending index
%     idx_end = idx - 1;
%     % get measurements data from current time step
%     X = dat.X(idx_start:idx_end);
%     Y = dat.Y(idx_start:idx_end);
%     Z = dat.Z(idx_start:idx_end);
%     B = dat.B(idx_start:idx_end);
%     prange = dat.rho(idx_start:idx_end);
%     % solve for position at current time step
%     sol = solve_pos(x0, X, Y, Z, B, prange);
%     x(:,tstep) = sol(1:3);
%     
% end

% solve with all pulled QOIs
while QOIs.Rx_time(idx) <= QOIs.Rx_time(end)-1 % set data pulling condition
    tstep = tstep + 1;
    % store current time step starting index
    idx_start = idx; 
    % check rows until time step changes
    time_inst = QOIs.Rx_time(idx_start);
    while QOIs.Rx_time(idx) == time_inst

        idx = idx+1;
    end
    % store current time step ending index
    idx_end = idx - 1;
    % get measurements data from current time step
    X = QOIs.X(idx_start:idx_end);
    Y = QOIs.Y(idx_start:idx_end);
    Z = QOIs.Z(idx_start:idx_end);
    B = QOIs.B(idx_start:idx_end);
    prange = QOIs.rho(idx_start:idx_end);
    sigma_rho = QOIs.sigma_rho(idx_start:idx_end);
    % calc WLS weight mtx
    W_rho = diag(1./sigma_rho);
    % solve for position at current time step
    sol = solve_pos(x0, X, Y, Z, B, prange,W_rho);
    x(:,tstep) = sol(1:3);
    
end

figure(10)
plot3(x(1,1:end),x(2,1:end),x(3,1:end),'-x');hold on
hold on
grid on
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)')
axis equal


%% functions

function rho0 = get_expected_pseudoranges(x_est, X, Y, Z, B)
rho0 = zeros(length(X),1);
for i=1:length(X)
    look_vec = [X(i);Y(i);Z(i)] - x_est(1:3);
    rho0(i) = norm(look_vec) + x_est(4) - B(i);
end

end


function G = get_geometry_matrix(x_est, X, Y, Z)
G = zeros(length(X), 4);
for i = 1:length(X)
    look_vec = [X(i);Y(i);Z(i)] - x_est(1:3);
    G(i,:) = [(-look_vec./norm(look_vec))' 1];
end

end

function sol = solve_pos(x0, X, Y, Z, B, prange,W_rho)
% X,Y,Z,B,prange col vecs wrt to each SV
x_est = x0;
res = 1e5;
it = 0;
while res > .05
    rho0 = get_expected_pseudoranges(x_est, X, Y, Z, B);
    drho = prange - rho0;
    G = get_geometry_matrix(x_est, X, Y, Z);
    dx = inv(G'*W_rho*G)*G'*W_rho*drho;
%     dx = inv(G'*G)*G'*drho;
    x_est = x_est + dx;
    res = norm(dx);
    it = it+1;
end

sol = x_est;
end


