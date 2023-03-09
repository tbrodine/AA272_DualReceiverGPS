
function [gnss,QOIs] = PreProcessGNSS(dirName,prFileName)
addpath('C:\MATLAB\Grad School\AA272\project/gpstools/opensource');

% download nasa ephem data to current folder

% prFileName = 'groundtruth.txt';
% dirName = 'C:\MATLAB\Grad School\AA272\project';
[ephem_struc, gnss, QOIs] = ProcessGnssMeasScript(dirName,prFileName);
ephem_raw = struct2table(ephem_struc);
gnss = struct2table(gnss);

ephem = renamevars(ephem_raw,["Asqrt", "Delta_n","e",           "OMEGA", "i0","OMEGA_DOT","GPS_Week","ttx"],...
                         ["sqrtA","DeltaN",  "Eccentricity","Omega0","Io","OmegaDot", "GPSWeek","TransTime"]);

%% Get X,Y,Z,B
% c = 2.99792458e8;
% 
% for idx = 1: length(gnss.Svid)
%     % calc SV ECEF for each measurement
%     SVID = gnss.Svid(idx);
%     eph_rows = ephem.PRN==SVID;
%     eph = ephem(eph_rows,:);
% %     [X_eph(idx,1),Y_eph(idx,1),Z_eph(idx,1)] = ...
% %         get_sat_ECEF(eph(1,:), gnss.ReceivedSvTimeNanos(idx)/1e9);
% %     B_eph(idx,1) = ...
% %         SV_clock_B(eph(1,:), gnss.ReceivedSvTimeNanos(idx)/1e9);
% %     idx = idx + 1;
% 
%     % or pull from post process
%     eph_raw = ephem_raw(eph_rows,:);
%     eph_struc = table2struct(eph_raw(1,:));
%     Nw = floor(-1.*double(gnss.FullBiasNanos(idx))/604800e9);
%     [xyzM,dtsvS] = GpsEph2Xyz(eph_struc,[Nw,gnss.ReceivedSvTimeNanos(idx)/1e9]);
%     X_eph(idx,1) = xyzM(1); Y_eph(idx,1) = xyzM(2); Z_eph(idx,1) = xyzM(3); 
%     B_eph(idx,1) = dtsvS*c;
% end
% gnss.RxTime_s = gnss.TimeNanos/1e9;
% gnss.X = X_eph; gnss.Y = Y_eph; gnss.Z = Z_eph; gnss.B = B_eph;

%% Get pseudorange

% Nw = floor(-1.*double(gnss.FullBiasNanos)./604800e9);
% Nw = 2251;
% t_Rx_hw = double(gnss.TimeNanos) + gnss.TimeOffsetNanos;
% b_hw = double(gnss.FullBiasNanos) + gnss.BiasNanos;
% t_Rx_GPS = t_Rx_hw - b_hw;
% t_Rx_w = t_Rx_GPS - Nw.*604800e9;
% rho_ns = t_Rx_w - double(gnss.ReceivedSvTimeNanos);
% rho = rho_ns.*299792458/1e9;
% gnss.rho = rho;


end

%% Functions

function [X,Y,Z] = get_sat_ECEF(ephem, tx_time)
mu = 3.986005e14; % from lecture 8
Omg_dot_E = 7.2921151467e-5; % from lecture 8

a = ephem.sqrtA^2;

n = sqrt(mu/a^3) + ephem.DeltaN;

t_k = tx_time - ephem.Toe;

M_k = double(ephem.M0 + n*t_k);

e = double(ephem.Eccentricity);
E_k = kepler(M_k,e);

sinvk = sqrt(1-e^2)*sin(E_k)/(1-e*cos(E_k));
cosvk = (cos(E_k) - e)/(1-e*cos(E_k));
v_k = atan2(sinvk,cosvk);

phi_k = v_k + ephem.omega;

dphi_k = ephem.Cus*sin(2*phi_k) + ephem.Cuc*cos(2*phi_k);
u_k = phi_k+dphi_k;

for l = 1:8
    dphi_k = ephem.Cus*sin(2*u_k) + ephem.Cuc*cos(2*u_k);
    u_k = phi_k+dphi_k; 
end


dr_k = ephem.Crs*sin(2*phi_k) + ephem.Crc*cos(2*phi_k);
di_k = ephem.Cis*sin(2*phi_k) + ephem.Cic*cos(2*phi_k);

Omg_k = double(ephem.Omega0 - Omg_dot_E*tx_time + ephem.OmegaDot*t_k);

r_k = a*(1-e*cos(E_k)) + dr_k;

i_k = double(ephem.Io + ephem.IDOT*t_k + di_k);

x_p = r_k*cos(u_k);
y_p = r_k*sin(u_k);

X = x_p*cos(Omg_k) - y_p*cos(i_k)*sin(Omg_k);
Y = x_p*sin(Omg_k) + y_p*cos(i_k)*cos(Omg_k);
Z = y_p*sin(i_k);

if ephem.IODE ~= ephem.IODC
    disp(['IODE err sat: ' num2str(ephem.PRN)])
end

end

function E = kepler(M,e)
dE = 1;
E = M;
while abs(dE) > 1e-6
    dE = -(E - e*sin(E)-M)/(1-e*cos(E));
    E = E + dE;
end
end

function B = SV_clock_B(ephem, tx_time)
% F = -4.442807633e-10;
c = 2.99792458e8;
mu = 3.986005e14;
F = -2*sqrt(mu)/c^2;

e = ephem.Eccentricity;
sqrtA = ephem.sqrtA;
T_0c = ephem.TransTime;
a_f0 = ephem.af0;
a_f1 = ephem.af1;
a_f2 = ephem.af2;
TGD = ephem.TGD;

% get eccentric anomaly E_k
a = sqrtA^2;
n = sqrt(mu/a^3) + ephem.DeltaN;
t_k = tx_time - ephem.Toe;
M_k = double(ephem.M0 + n*t_k);
e = ephem.Eccentricity;
E_k = kepler(M_k,e);

% relativistic term
dt_r = F*e*sqrtA*sin(E_k);

B = a_f0 + a_f1*double(tx_time-T_0c) + a_f2*double(tx_time-T_0c)^2 + dt_r - TGD;
B = B*c;

end