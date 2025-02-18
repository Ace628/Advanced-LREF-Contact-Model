clc
close all
clear all

%% Material properties
%material properties of laminates
%mat(id#,:)=[E11 E22 E33 v12 v13 v23 G12 G13 G23]; % MPa
% mat(1,:)=[20391.63 28.52 28.52 0.39 0.39 0.52 9.64 9.64 9.37];   % Fiber-reinforced rubber
% mat(1,:)=[20384.86 13.06 13.06 0.389 0.389 0.898 3.50 3.50 3.43]; % Fiber-reinforced rubber
% mat(2,:)=[15 15 15 0.48 0.48 0.48 5.068 5.068 5.068]; % RSM rubber
mat(1,:)=[20000 20 20 0.4 0.4 0.6667 10 10 6];
E_ex = 4.8;
nu_ex = 0.48;

%% Composites
n = 5;         % Ply layers
R = 200;       % The radius of the middle surface of the ring
h = 20;        % The thickness of the whole ring [mm]
b = 60;        % The width of the ring [mm]
Plyt = h/n;    % Ply thickness [mm]
R_in = R-h/2;
Neutral_R = R;
he = 8;
R_ext = R_in +h;

%Layup properties [Material ID, ply thickness, orientation (deg)]
Layup(1,:)=[1,Plyt,90];
Layup(2,:)=[1,Plyt,0];
Layup(3,:)=[1,Plyt,90];
Layup(4,:)=[1,Plyt,0];
Layup(5,:)=[1,Plyt,90];

% Layup(1,:)=[1,Plyt,0];
% Layup(2,:)=[1,Plyt,90];
% Layup(3,:)=[1,Plyt,90];
% Layup(4,:)=[1,Plyt,0];


%% Neutral axis & ABD matrix
[threeDQbar] = Qbar(Layup,mat);
[z]=mid_z(Layup,h);
[A,B,D,A55]=ABD_curved(z,Neutral_R,b,h,threeDQbar,Layup);
[L_ex, ABD_ex]=ABD_curved_ex(Neutral_R,b,h,he,E_ex,nu_ex);

%% Define the foundation and contact stiffnesses
K_R = 1;       % The radial stiffness of the elastic foundation [N/mm]

% Unilateral foundation case
K_R_C = 1e-1*K_R;  % The compressive stiffness of the elastic foundation [N/mm]
K_T = 1e-1*K_R;% The tangential stiffness of the elastic foundation [N/mm]

k_gen = 0.2; % The contact stiffness

%% Nondimensional parameters related to material and geometric parameters
[U_Cons,Roots,Delta,UTh_Cons,UFi_Cons,UPsi_Cons,K_ring] = Solving_Constants(A,B,D,A55,L_ex,ABD_ex,Neutral_R,h,K_R,K_T);

%% The stiffness of the ring-foundation system
Delta_initial = ur_bar_Laminated_REF(0,0,U_Cons,Roots,Delta);
K_LRLF = 1/Delta_initial;
Nor_K_LRLF = K_LRLF/K_ring;

Point_F = 5*K_ring;
D_ang_PointF = pi/1800;
ang_PointF0 = [pi/2:D_ang_PointF:pi,-pi+D_ang_PointF:D_ang_PointF:-pi/2]';

[q_Cps_PointF,Rot_Cps_PointF,K_LRNF] = Compensation_Load_NF_PointLoad3(R,h,K_R,K_R_C,Point_F,ang_PointF0,D_ang_PointF,U_Cons,Roots,Delta,UTh_Cons,UFi_Cons);
Nor_K_LRNF = K_LRNF/K_ring;

%% Contact Behaviors of the ring on a nonreciprocal foundation
delta0 = 20;  % Vertical deflection of the rigid surface

% Plate contact
N_Cnt = 400;
Cnt_ang0 = acos(1-delta0/R_ext);
D_Cnt_ang = 2*Cnt_ang0/N_Cnt;
ang_rot = (-Cnt_ang0:D_Cnt_ang:Cnt_ang0)';              % The rotation angle within the penetration region 
ang_i = [pi-Cnt_ang0:D_Cnt_ang:pi,-pi+D_Cnt_ang:D_Cnt_ang:-pi+Cnt_ang0]'; % The initial penetration region

% Defining variables used in the contact behavior
Cnt_ur0 = zeros(length(ang_i),1);
Cnt_ut0 = zeros(length(ang_i),1);
Cnt_fi0 = zeros(length(ang_i),1);
Cnt_ur = zeros(length(ang_i),1);
Cnt_ut = zeros(length(ang_i),1);
Cnt_fi = zeros(length(ang_i),1);
Cnt_UR = zeros(length(ang_i),1);

Cnt_urF = zeros(length(ang_i),1);
Cnt_utF = zeros(length(ang_i),1);
Cnt_fiF = zeros(length(ang_i),1);

D_Fr = zeros(length(ang_i),1);
D_Fr0 = zeros(length(ang_i),1);
F_i = zeros(length(ang_i),1);
Distance = zeros(length(ang_i),1);
r_deform = R_ext*ones(length(ang_i),1);
theta_di = ang_i;
rot_cnt = ang_i;

Road = zeros(length(ang_rot),1);

% Defining variables of the whole ring used for compendation load
D_ang_Cps = pi/1000;
% ang_Cps0 = (-pi:D_ang_Cps:pi)';
ang_Cps0 = [pi/2:D_ang_Cps:pi,-pi+D_ang_Cps:D_ang_Cps:-pi/2]';

Cnt_urCps = zeros(length(ang_Cps0),1);
Cnt_utCps = zeros(length(ang_Cps0),1);
Cnt_fiCps = zeros(length(ang_Cps0),1);

% Contact behaviors in the contact region
for j0 = 1:20000 % The maximum iteration step

for i1 = 1:length(ang_i)
    Road(i1) = (R+h/2-delta0)/cos(ang_rot(i1));   % Road function for calculating geometry errors
    Distance(i1) = r_deform(i1) - Road(i1);
    D_Fr0(i1) = k_gen*Distance(i1);               % The distributing force along the contact region
    Delta_F = F_i(i1) + D_Fr0(i1);
    if Delta_F > 0
        D_Fr(i1) = D_Fr0(i1);                     % Ensuring the direction of the pressure
    else
        D_Fr(i1) = 0;
    end
end

F_i = F_i + D_Fr;                                 % Collecting the total distributed load p_c=F_i/b

% The compensation load under F_i in each iterative step
[q_Cps,Rot_Cps] = Compensation_Load_NF2(R,h,K_R,K_R_C,D_Cnt_ang,F_i,ang_rot,ang_Cps0,D_ang_Cps,U_Cons,Roots,Delta,UTh_Cons,UFi_Cons);

% The total deformation of the penetration region under F_i and q_Cps
for i2 = 1:length(ang_i)

    for i3 = 1:length(F_i)
        Cnt_urF(i3) = F_i(i3)*(R+h/2)*D_Cnt_ang*ur_bar_Laminated_REF(ang_i(i2),ang_rot(i3),U_Cons,Roots,Delta);
        Cnt_utF(i3) = F_i(i3)*(R+h/2)*D_Cnt_ang*ut_bar_Laminated_REF(ang_i(i2),ang_rot(i3),U_Cons,Roots,Delta,UTh_Cons);
        Cnt_fiF(i3) = F_i(i3)*(R+h/2)*D_Cnt_ang*ufi_bar_Laminated_REF(ang_i(i2),ang_rot(i3),U_Cons,Roots,Delta,UFi_Cons);
    end

    for i4 = 1:length(q_Cps)
        Cnt_urCps(i4) = q_Cps(i4)*(R+h/2)*D_ang_Cps*ur_bar_Laminated_REF(ang_i(i2),Rot_Cps(i4),U_Cons,Roots,Delta);
        Cnt_utCps(i4) = q_Cps(i4)*(R+h/2)*D_ang_Cps*ut_bar_Laminated_REF(ang_i(i2),Rot_Cps(i4),U_Cons,Roots,Delta,UTh_Cons);
        Cnt_fiCps(i4) = q_Cps(i4)*(R+h/2)*D_ang_Cps*ufi_bar_Laminated_REF(ang_i(i2),Rot_Cps(i4),U_Cons,Roots,Delta,UFi_Cons);
    end

    Cnt_ur(i2) = sum(Cnt_urF)+sum(Cnt_urCps);
    Cnt_ut(i2) = sum(Cnt_utF)+sum(Cnt_utCps);
    Cnt_fi(i2) = sum(Cnt_fiF)+sum(Cnt_fiCps);

    r_deform(i2) = sqrt((R+h/2+Cnt_ur(i2))^2+(Cnt_ut(i2)+h/2*Cnt_fi(i2))^2);
    theta_di(i2) = ang_i(i2) + atan((Cnt_ut(i2)+h/2*Cnt_fi(i2))/(R+h/2+Cnt_ur(i2)));
    rot_cnt(i2) = theta_di(i2)-sign(theta_di(i2))*pi;
end

% The error between the road and the penetration region
D_err = Cnt_ur - Cnt_UR;

% Stop the iteration if the change of the deformation is quite small
if max(abs(D_err)) < 4e-6      %4e-7 2.5e-6 Predefined threshold
    break
else
    Cnt_UR = Cnt_ur;

    Cnt_ur = zeros(length(ang_i),1);
    Cnt_ut = zeros(length(ang_i),1);
    Cnt_fi = zeros(length(ang_i),1);
    continue
end

end

% The vertical reaction force on the plate
F_v = zeros(length(F_i),1);

for i6 = 1:length(F_i)
    F_v(i6) = F_i(i6)*cos(theta_di(i6))*(R+h/2)*D_Cnt_ang;
end

F_ver = sum(F_v);

%% Plotting the results
D_ang_j = pi/1800;
ang_j = (-pi:D_ang_j:pi)';

Cnt_UR1 = zeros(length(ang_j),1);
Cnt_UT1 = zeros(length(ang_j),1);
Cnt_FI1 = zeros(length(ang_j),1);

Cnt_UR2 = zeros(length(F_i),1);
Cnt_UT2 = zeros(length(F_i),1);
Cnt_FI2 = zeros(length(F_i),1);

Cnt_UR3 = zeros(length(q_Cps),1);
Cnt_UT3 = zeros(length(q_Cps),1);
Cnt_FI3 = zeros(length(q_Cps),1);

R_deform = zeros(length(ang_j),1);
Theta_Di = zeros(length(ang_j),1);

F_kv0 = zeros(length(ang_j),1);

R1 = R_ext*ones(length(ang_j),1);

for j1 = 1:length(ang_j)
    
    for j2 = 1:length(F_i)
        Cnt_UR2(j2) = F_i(j2)*R_ext*D_Cnt_ang*ur_bar_Laminated_REF(ang_j(j1),ang_rot(j2),U_Cons,Roots,Delta);
        Cnt_UT2(j2) = F_i(j2)*R_ext*D_Cnt_ang*ut_bar_Laminated_REF(ang_j(j1),ang_rot(j2),U_Cons,Roots,Delta,UTh_Cons);
        Cnt_FI2(j2) = F_i(j2)*R_ext*D_Cnt_ang*ufi_bar_Laminated_REF(ang_j(j1),ang_rot(j2),U_Cons,Roots,Delta,UFi_Cons);
    end

    for j3 = 1:length(q_Cps)
        Cnt_UR3(j3) = q_Cps(j3)*R_ext*D_ang_Cps*ur_bar_Laminated_REF(ang_j(j1),Rot_Cps(j3),U_Cons,Roots,Delta);
        Cnt_UT3(j3) = q_Cps(j3)*R_ext*D_ang_Cps*ut_bar_Laminated_REF(ang_j(j1),Rot_Cps(j3),U_Cons,Roots,Delta,UTh_Cons);
        Cnt_FI3(j3) = q_Cps(j3)*R_ext*D_ang_Cps*ufi_bar_Laminated_REF(ang_j(j1),Rot_Cps(j3),U_Cons,Roots,Delta,UFi_Cons);
    end

    Cnt_UR1(j1) = sum(Cnt_UR2) + sum(Cnt_UR3);
    Cnt_UT1(j1) = sum(Cnt_UT2) + sum(Cnt_UT3);
    Cnt_FI1(j1) = sum(Cnt_FI2) + sum(Cnt_FI3);

    R_deform(j1) = sqrt((R+h/2+Cnt_UR1(j1))^2+(Cnt_UT1(j1)+h/2*Cnt_FI1(j1))^2);
    Theta_Di(j1) = ang_j(j1) + atan((Cnt_UT1(j1)+h/2*Cnt_FI1(j1))/(R+h/2+Cnt_UR1(j1)));

    if Cnt_UR1(j1) < 0
        F_kv0(j1) = (K_T*(R-h/2)*(Cnt_UT1(j1)-h/2*Cnt_FI1(j1))*sin(ang_j(j1)))*D_ang_j;
    else
        F_kv0(j1) = (-K_R*(R-h/2)*Cnt_UR1(j1)*cos(ang_j(j1))+K_T*(R-h/2)*(Cnt_UT1(j1)-h/2*Cnt_FI1(j1))*sin(ang_j(j1)))*D_ang_j;
    end
end

F_kv = sum(F_kv0);     % Calculation of the vertical reaction force on the elastic foundation

figure(1)
polarplot(Theta_Di,R_deform,'r')  % Deformed ring due to the distributed load
hold on

% polarplot(ang_j,R1,'k')        % Undeformed ring
polarplot(ang_i,Road)

pax = gca;
pax.ThetaDir = 'clockwise';
pax.ThetaZeroLocation = 'top';

hold off

figure(2)
plot(rad2deg(rot_cnt),F_i/b)

% plot((R+h/2)*sin(ang_i),F_i)

% % legend('Inextensible','Extensible','Undeformed')
% % legend('Linear','Nonlinear','Undeformed')
