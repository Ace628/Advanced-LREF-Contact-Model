clc
close all
clear all

%% Material properties
% Material properties of fibers and matrix
E_f = 40000;   % Elastic modulus of fibers [MPa]
nu_f = 0.22;   % Poisson's ratio of fibers
E_m = 5.0;    % Elastic modulus of matrix [MPa] eta=1 E_m=4.834
nu_m = 0.48;   % Poisson's ratio of matrix
V_f = 0.5;     % The volume fraction of fibers

%material properties of laminates
%mat(id#,:)=[E11 E22 E33 v12 v13 v23 G12 G13 G23]; % MPa
% mat(1,:)=[20391.63 28.52 28.52 0.39 0.39 0.52 9.64 9.64 9.37]; % Fiber-reinforced rubber
% mat(1,:)=[20384.86 13.06 13.06 0.389 0.389 0.898 3.50 3.50 3.43];
% mat(1,:)=Self_Consistant_Field(E_f,nu_f,E_m,nu_m,V_f);
% mat(1,:)=[6740 1142.61 1142.61 0.45 0.45 0.48 386.23 386.23 385.80];
% mat(1,:)=[20000 20 20 0.4 0.4 0.6667 10 10 6];
% mat(1,:)=[20000 20000 20000 0.4 0.4 0.4 7142.857 7142.857 7142.857];
% mat(2,:)=[15 15 15 0.48 0.48 0.48 5.068 5.068 5.068]; % RSM rubber
E_ex = 5;
nu_ex = 0.48;

% G/E cases
% mat(1,:) = [22908.46 6701.13 6701.13 0.39 0.39 0.52 2290.87 2290.87 2200.89]; %k=3.060 G/E=1e-1
% mat(1,:) = [21532.5 3168.72 3168.72 0.39 0.39 0.52 1076.63 1076.63 1040.72];  %k=2.450 G/E=5e-2
% mat(1,:) = [20596.95 609.03 609.03 0.39 0.39 0.52 205.98 205.98 200.03];      %k=1.965 G/E=1e-2
% mat(1,:) = [20488.43 303.09 303.09 0.39 0.39 0.52 102.45 102.45 99.54];       %k=1.800 G/E=5e-3
% mat(1,:) = [20402.85 60.44 60.44 0.39 0.39 0.52 20.42 20.42 19.85];           %k=1.155 G/E=1e-3
% mat(1,:) = [20392.21 30.19 30.19 0.39 0.39 0.52 10.2 10.2 9.92];              %k=0.800 G/E=5e-4
% mat(1,:) = [20383.72 6.04 6.04 0.39 0.39 0.52 2.04 2.04 1.98];                %k=0.233 G/E=1e-4
% mat(1,:) = [20382.66 3.02 3.02 0.39 0.39 0.52 1.02 1.02 0.99];                %k=0.126 G/E=5e-5
mat(1,:) = [20381.81 0.6 0.6 0.39 0.39 0.52 0.2 0.2 0.2];                     %k=0.0318 G/E=1e-5

%% Composites
n = 5;    % Ply layers
R = 372;    % The radius of the middle surface of the ring
h = 20;      % The thickness of the whole ring [mm]
R_in = R-h/2;
b = 168;       % The width of the ring [mm]
Plyt = h/n;    % Ply thickness [mm]
Neutral_R = R;
R_ext = R+h/2;
he = 8;

%Layup properties [Material ID, ply thickness, orientation (deg)]
% Layup(1,:)=[1,Plyt,90];
% Layup(2,:)=[1,Plyt,0];
% Layup(3,:)=[1,Plyt,90];
% Layup(4,:)=[1,Plyt,0];
% Layup(5,:)=[1,Plyt,90];

Layup(1,:)=[1,Plyt,0];
Layup(2,:)=[1,Plyt,90];
Layup(3,:)=[1,Plyt,90];
Layup(4,:)=[1,Plyt,90];
Layup(5,:)=[1,Plyt,0];

%% ABD matrix
[threeDQbar] = Qbar(Layup,mat);
[z]=mid_z(Layup,h);
[A,B,D,A55]=ABD_curved(z,Neutral_R,b,h,threeDQbar,Layup);
[L_ex, ABD_ex]=ABD_curved_ex(Neutral_R,b,h,he,E_ex,nu_ex);

%% effective stiffnesses and nondimensional parameters
K_R = 0.0318;        % The radial foundation stiffness [N/mm]
K_T = 1e-3*K_R;      % The tangential foundatin stiffness [N/mm]
k_gen = 0.1*K_R;     % The assumed contact stiffness

[U_Cons,Roots,Delta,UTh_Cons,UFi_Cons,UPsi_Cons,K_ring] = Solving_Constants(A,B,D,A55,L_ex,ABD_ex,Neutral_R,h,K_R,K_T);

%% The stiffness of the ring-foundation system
Delta_initial = ur_bar_Laminated_REF(0,0,U_Cons,Roots,Delta);
K_LRLF = 1/Delta_initial;
Nor_K_LRLF = K_LRLF/K_ring;

%% Contact responses
% Plate contact
delta0 = 0.08*R;                      % Vertical deflection of the rigid surface
N_Cnt = 400;
Cnt_ang0 = acos(1-delta0/R_ext);
D_Cnt_ang = 2*Cnt_ang0/N_Cnt;
ang_rot = (-Cnt_ang0:D_Cnt_ang:Cnt_ang0)';              % The rotation angle within the penetration region 
theta_di0 = [pi-Cnt_ang0:D_Cnt_ang:pi,-pi+D_Cnt_ang:D_Cnt_ang:-pi+Cnt_ang0]'; % The initial penetration region
r_deform0 = (R+h/2)*ones(length(theta_di0),1);
F_i = zeros(length(theta_di0),1);

[F_i_C,rot_cnt,F_ver,j0,r_deform_C,theta_di_C] = PlateContact(k_gen,R,h,delta0,D_Cnt_ang,r_deform0,theta_di0,ang_rot,U_Cons,Roots,Delta,UTh_Cons,UFi_Cons);

%% Plotting the results  
 Road0 = (R+h/2-delta0)./cos(ang_rot);

D_ang_j = pi/1800;
ang_j = (-pi:D_ang_j:pi)';

Cnt_UR1 = zeros(length(ang_j),1);
Cnt_UT1 = zeros(length(ang_j),1);
Cnt_FI1 = zeros(length(ang_j),1);
Cnt_UR2 = zeros(length(F_i_C),1);
Cnt_UT2 = zeros(length(F_i_C),1);
Cnt_FI2 = zeros(length(F_i_C),1);

R_deform = zeros(length(ang_j),1);
Theta_Di = zeros(length(ang_j),1);

F_kv0 = zeros(length(ang_j),1);

R1 = R_ext*ones(length(ang_j),1);

for j1 = 1:length(ang_j)
    
    for j2 = 1:length(F_i_C)
        Cnt_UR2(j2) = F_i_C(j2)*R_ext*D_Cnt_ang*ur_bar_Laminated_REF(ang_j(j1),rot_cnt(j2),U_Cons,Roots,Delta);
        Cnt_UT2(j2) = F_i_C(j2)*R_ext*D_Cnt_ang*ut_bar_Laminated_REF(ang_j(j1),rot_cnt(j2),U_Cons,Roots,Delta,UTh_Cons);
        Cnt_FI2(j2) = F_i_C(j2)*R_ext*D_Cnt_ang*ufi_bar_Laminated_REF(ang_j(j1),rot_cnt(j2),U_Cons,Roots,Delta,UFi_Cons);
    end

    Cnt_UR1(j1) = sum(Cnt_UR2);
    Cnt_UT1(j1) = sum(Cnt_UT2);
    Cnt_FI1(j1) = sum(Cnt_FI2);

    R_deform(j1) = sqrt((R+h/2+Cnt_UR1(j1))^2+(Cnt_UT1(j1)+h/2*Cnt_FI1(j1))^2);
    Theta_Di(j1) = ang_j(j1) + atan((Cnt_UT1(j1)+h/2*Cnt_FI1(j1))/(R+h/2+Cnt_UR1(j1)));

    F_kv0(j1) = (-K_R*Cnt_UR1(j1)*cos(Theta_Di(j1))+K_T*(Cnt_UT1(j1)-h/2*Cnt_FI1(j1))*sin(Theta_Di(j1)))*(R-h/2)*D_ang_j;
end

F_kv = sum(F_kv0);     % Calculation of the vertical reaction force on the elastic foundation

figure(1)
polarplot(Theta_Di,R_deform,'r')  % Deformed ring due to the distributed load
pax = gca;
pax.ThetaDir = 'clockwise';
pax.ThetaZeroLocation = 'top';
hold on
polarplot(ang_j,R1,'--k')         % Undeformed ring
polarplot(theta_di0,Road0,'--k')   % Road
hold off

figure(2)
plot(rad2deg(rot_cnt),F_i_C/b);

% plot((R+h/2)*sin(ang_i),F_i)

% % legend('Inextensible','Extensible','Undeformed')
% % legend('Linear','Nonlinear','Undeformed')
