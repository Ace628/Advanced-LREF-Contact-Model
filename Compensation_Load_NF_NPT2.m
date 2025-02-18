function [q_Cps,Rot_Cps] = Compensation_Load_NF_NPT2(R,h,he,K_R,K_R_C,D_Cnt_ang,F_i,ang_rot,ang_j,D_ang_j,U_Cons,Roots,Delta,UTh_Cons,UFi_Cons,UPsi_Cons)

D_ang_Cps = D_ang_j;

Cps_ur0 = zeros(length(F_i),1);
Cps_ut0 = zeros(length(F_i),1);
Cps_fi0 = zeros(length(F_i),1);
Cps_psi0 = zeros(length(F_i),1);

Cps_ur1 = zeros(length(ang_j),1);
Cps_ut1 = zeros(length(ang_j),1);
Cps_fi1 = zeros(length(ang_j),1);
Cps_psi1 = zeros(length(ang_j),1);

Theta_Di_Cps0 = zeros(length(ang_j),1);
R_deform_Cps0 = zeros(length(ang_j),1);

D_q_Cps0 = zeros(length(ang_j),1);
Rot_Cps0 = zeros(length(ang_j),1);

% Getting the initial compensation load
for i7 = 1:length(ang_j)

    for i8 = 1:length(F_i)
        Cps_ur0(i8) = F_i(i8)*(R+h/2+he)*D_Cnt_ang*ur_bar_Laminated_REF(ang_j(i7),ang_rot(i8),U_Cons,Roots,Delta);
        Cps_ut0(i8) = F_i(i8)*(R+h/2+he)*D_Cnt_ang*ut_bar_Laminated_REF(ang_j(i7),ang_rot(i8),U_Cons,Roots,Delta,UTh_Cons);
        Cps_fi0(i8) = F_i(i8)*(R+h/2+he)*D_Cnt_ang*ufi_bar_Laminated_REF(ang_j(i7),ang_rot(i8),U_Cons,Roots,Delta,UFi_Cons);
        Cps_psi0(i8) = F_i(i8)*(R+h/2+he)*D_Cnt_ang*upsi_bar_Laminated_REF(ang_j(i7),ang_rot(i8),U_Cons,Roots,Delta,UPsi_Cons);
    end
    
    Cps_ur1(i7) = sum(Cps_ur0);
    Cps_ut1(i7) = sum(Cps_ut0);
    Cps_fi1(i7) = sum(Cps_fi0);
    Cps_psi1(i7) = sum(Cps_psi0);

    Theta_Di_Cps0(i7) = ang_j(i7) + atan((Cps_ut1(i7)+(h/2+he)*Cps_fi1(i7))/((R+h/2+he)+Cps_ur1(i7)+he*Cps_psi1(i7)));
    R_deform_Cps0(i7) = sqrt((R+h/2+he+Cps_ur1(i7)+he*Cps_psi1(i7))^2+(Cps_ut1(i7)+(he+h/2)*Cps_fi1(i7))^2);

    if Cps_ur1(i7) < 0
        D_q_Cps0(i7) = -(K_R-K_R_C)*Cps_ur1(i7);
    else
        D_q_Cps0(i7) = 0;
    end
    
    Rot_Cps0(i7) = Theta_Di_Cps0(i7) - sign(Theta_Di_Cps0(i7))*pi;        
end


Cps_Ur = Cps_ur1;
Cps_Ut = Cps_ut1;
Cps_Fi = Cps_fi1;
Cps_Psi = Cps_psi1;

D_q_Cps = D_q_Cps0;
Rot_Cps = Rot_Cps0;

Theta_Di_Cps = Theta_Di_Cps0;
% R_deform_Cps = R_deform_Cps0;

q_Cps = zeros(length(Theta_Di_Cps0),1);

Cps_ur = zeros(length(q_Cps),1);
Cps_ut = zeros(length(q_Cps),1);
Cps_fi = zeros(length(q_Cps),1);
Cps_psi = zeros(length(q_Cps),1);
Cps_ur2 = zeros(length(q_Cps),1);
Cps_ut2 = zeros(length(q_Cps),1);
Cps_fi2 = zeros(length(q_Cps),1);
Cps_psi2 = zeros(length(q_Cps),1);

% R_deform_Cps2 = R_deform_Cps0;

% Getting the total deformation
for c0 = 1:50000
    
    D_q_Cps1 = D_q_Cps;
    q_Cps = q_Cps + D_q_Cps;
    D_q_Cps = zeros(length(q_Cps),1);
%     Theta_Di_Cps1 = zeros(length(Theta_Di_Cps0),1);
%     R_deform_Cps1 = zeros(length(Theta_Di_Cps0),1);

% Deformation under each increased compensation load
    for c1 = 1:length(Theta_Di_Cps)

        for c2 = 1:length(D_q_Cps1)
            Cps_ur2(c2) = D_q_Cps1(c2)*D_ang_Cps*(R+h/2+he)*ur_bar_Laminated_REF(Theta_Di_Cps(c1),Rot_Cps(c2),U_Cons,Roots,Delta);
            Cps_ut2(c2) = D_q_Cps1(c2)*D_ang_Cps*(R+h/2+he)*ut_bar_Laminated_REF(Theta_Di_Cps(c1),Rot_Cps(c2),U_Cons,Roots,Delta,UTh_Cons);
            Cps_fi2(c2) = D_q_Cps1(c2)*D_ang_Cps*(R+h/2+he)*ufi_bar_Laminated_REF(Theta_Di_Cps(c1),Rot_Cps(c2),U_Cons,Roots,Delta,UFi_Cons);
            Cps_psi2(c2) = D_q_Cps1(c2)*D_ang_Cps*(R+h/2+he)*upsi_bar_Laminated_REF(Theta_Di_Cps(c1),Rot_Cps(c2),U_Cons,Roots,Delta,UPsi_Cons);
        end

        Cps_ur(c1) = sum(Cps_ur2);
        Cps_ut(c1) = sum(Cps_ut2);
        Cps_fi(c1) = sum(Cps_fi2);
        Cps_psi(c1) = sum(Cps_psi2);

        Cps_Ur(c1) = Cps_Ur(c1) + Cps_ur(c1);
        Cps_Ut(c1) = Cps_Ut(c1) + Cps_ut(c1);
        Cps_Fi(c1) = Cps_Fi(c1) + Cps_fi(c1);
        Cps_Psi(c1) = Cps_Psi(c1) + Cps_psi(c1);

%         Theta_Di_Cps1(c1) = Theta_Di_Cps(c1) + atan((Cps_Ut(c1)+(h/2+he)*Cps_Fi(c1))/(R_deform_Cps(c1)+Cps_Ur(c1)+he*Cps_Psi(c1)));
%         R_deform_Cps1(c1) = sqrt((R_deform_Cps(c1)+Cps_Ur(c1)+he*Cps_Psi(c1))^2+(Cps_Ut(c1)+(h/2+he)*Cps_Fi(c1))^2); 

        if Cps_ur(c1) < 0
           D_q_Cps(c1) = -(K_R-K_R_C)*Cps_ur(c1);
        else
           D_q_Cps(c1) = 0;
        end

    end

% Stop the iteration if the change of the deformation is quite small
    diff = max(abs(Cps_ur));
    if diff < 8e-7
       break
    else
%         R_deform_Cps2 = R_deform_Cps1;

        Cps_ur = zeros(length(Theta_Di_Cps),1);
        Cps_ut = zeros(length(Theta_Di_Cps),1);
        Cps_fi = zeros(length(Theta_Di_Cps),1);
        Cps_psi = zeros(length(Theta_Di_Cps),1);

        continue
    end

end

end