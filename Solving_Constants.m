function [U_Cons,Roots,Delta,UTh_Cons,UFi_Cons,UPsi_Cons,K_ring] = Solving_Constants(A,B,D,A55,L_ex,ABD_ex,Neutral_R,h,K_R,K_T)

A11 = A(1,1);
B11 = B(1,1);
D11 = D(1,1);
Cab1 = ABD_ex(1) + ABD_ex(2)/(Neutral_R+h/2);
Cbd1 = ABD_ex(2) + ABD_ex(3)/(Neutral_R+h/2);
Cl12 = L_ex(1)*(Neutral_R+h/2) + L_ex(2);
Ckr = K_R*(Neutral_R-h/2);
Ckt = K_T*(Neutral_R-h/2);

EA = A(1,1)-B(1,1)^2/D(1,1);
EV = A(1,1)*D(1,1)/B(1,1)-B(1,1);
EI = D(1,1)-B(1,1)^2/A(1,1);
GA = A55;

K_ring = 8*pi/(pi^2*Neutral_R/GA+pi^2*(Neutral_R/EA+Neutral_R^2/EV)+(pi^2-8-8*EI/Neutral_R/EV)*(Neutral_R^3/EI+Neutral_R^2/EV));

C_r1 = 2-Neutral_R*Ckr/A55-Neutral_R*Ckt*(h^2*A11/4+h*B11+D11)/(A11*D11-B11^2);
C_r2 = 1+Neutral_R*Ckr*(Neutral_R^2*A11+2*Neutral_R*B11+D11)/(A11*D11-B11^2)+ ...
    Neutral_R*Ckt*(Neutral_R-h/2)*(h*A11+2*B11)/(A11*D11-B11^2)+ ...
    Neutral_R*Ckt/A55*(1+Neutral_R*Ckr*(h^2*A11/4+h*B11+D11)/(A11*D11-B11^2));
C_r3 = -Neutral_R*Ckt*(Neutral_R-h/2)^2*(Neutral_R*Ckr+A11)/(A11*D11-B11^2);

P = A55/(B11-Neutral_R*A55)*(1/A55+(Neutral_R^2*A11+2*Neutral_R*B11+D11)/(A11*D11-B11^2))+...
    Neutral_R*Ckt/(A11*D11-B11^2)*(h/2+h^2*(A11+A55)/4/(B11-Neutral_R*A55))+...
    Neutral_R*Ckt/(Neutral_R*A11+B11)*(1/A55+(Neutral_R*B11+D11)/(A11*D11-B11^2))*(1+h*(A11+A55)/2/(B11-Neutral_R*A55));

C_theta1 = 1/P/(Neutral_R*A11+B11);
C_theta2 = ((B11*P+1)*A55-Ckr*Neutral_R)/P/A55/(Neutral_R*A11+B11)-...
    Neutral_R*Ckt*(h*A11/2+B11)/(A11*D11-B11^2)/(Neutral_R*A11+B11)*(h/2+(A11*D11-B11^2)/A55/(Neutral_R*A11+B11)+(Neutral_R*B11+D11)/(Neutral_R*A11+B11));
C_theta3 = -(Neutral_R*Ckr+A11)/P/(Neutral_R*A11+B11)* ...
(Neutral_R*A55*(Neutral_R^2*A11+2*Neutral_R*B11+D11)/(A11*D11-B11^2)/(B11-Neutral_R*A55)+ ...
    (h*Neutral_R*Ckt/2+Neutral_R*A55)/A55/(B11-Neutral_R*A55)+ ...
    (h^2*Neutral_R*Ckt/4*(Neutral_R*A11+B11)+h*Neutral_R*Ckt/2*(Neutral_R*B11+D11))/(A11*D11-B11^2)/(B11-Neutral_R*A55));

C_phi1 = -((A11*D11-B11^2)+A55*(Neutral_R*B11+D11))/(Neutral_R*A55*(Neutral_R*B11+D11)+h^2*Neutral_R*Ckt*B11/4+h*Neutral_R*Ckt*D11/2);
C_phi2 = -(A11*D11-B11^2)/(Neutral_R*A55*(Neutral_R*B11+D11)+h^2*Neutral_R*Ckt*B11/4+h*Neutral_R*Ckt*D11/2);
C_phi3 = (A55*(Neutral_R*B11+D11)+h*Neutral_R*Ckt*B11/2+Neutral_R*Ckt*D11)/(Neutral_R*A55*(Neutral_R*B11+D11)+h^2*Neutral_R*Ckt*B11/4+h*Neutral_R*Ckt*D11/2);

C_psi1 = -Cab1/(Cl12+Cbd1);
C_psi2 = -Cbd1/(Cl12+Cbd1);

DA = C_r1^2-3*C_r2;
DB = C_r1*C_r2-9*C_r3;
DC = C_r2^2-3*C_r1*C_r3;

Delta = DB^2-4*DA*DC;

if Delta<0

    Q0 = (2*DA*C_r1-3*DB)/2/sqrt(DA^3);
    eta = acos(Q0);
    alpha1 = sqrt((-C_r1+sqrt(DA)*(cos(eta/3)+sqrt(3)*sin(eta/3)))/3);
    beta1 = sqrt((-C_r1+sqrt(DA)*(cos(eta/3)-sqrt(3)*sin(eta/3)))/3);
    gamma1 = sqrt((-C_r1-2*sqrt(DA)*cos(eta/3))/3);

    C_utheta1 = alpha1^3*C_theta1+alpha1*C_theta2+C_theta3/alpha1;
    C_utheta2 = beta1^3*C_theta1+beta1*C_theta2+C_theta3/beta1;
    C_utheta3 = gamma1^3*C_theta1+gamma1*C_theta2+C_theta3/gamma1;
    C_uphi1 = alpha1*C_phi1+alpha1^2*C_utheta1*C_phi2+C_utheta1*C_phi3;
    C_uphi2 = beta1*C_phi1+beta1^2*C_utheta2*C_phi2+C_utheta2*C_phi3;
    C_uphi3 = gamma1*C_phi1+gamma1^2*C_utheta3*C_phi2+C_utheta3*C_phi3;
    C_upsi1 = C_psi1+alpha1*C_utheta1*C_psi1+alpha1*C_uphi1*C_psi2;
    C_upsi2 = C_psi1+beta1*C_utheta2*C_psi1+beta1*C_uphi2*C_psi2;
    C_upsi3 = C_psi1+gamma1*C_utheta3*C_psi1+gamma1*C_uphi3*C_psi2;

    D0 = [alpha1 beta1 gamma1;C_utheta1 C_utheta2 C_utheta3;C_uphi1 C_uphi2 C_uphi3];
    HD1 = [C_utheta2 C_utheta3;C_uphi2 C_uphi3];
    HD2 = [C_utheta1 C_utheta3;C_uphi1 C_uphi3];
    HD3 = [C_utheta1 C_utheta2;C_uphi1 C_uphi2];
    H1 = det(HD1)/det(D0);
    H2 = det(HD2)/det(D0);
    H3 = det(HD3)/det(D0);

    C1 = -Neutral_R*H1/2/A55/sinh(alpha1*pi);
    C3 = Neutral_R*H2/2/A55/sinh(beta1*pi);
    C5 = -Neutral_R*H3/2/A55/sinh(gamma1*pi);

    U_Cons = [C1;C3;C5];
    Roots = [alpha1;beta1;gamma1];

    UTh_Cons = [C_utheta1;C_utheta2;C_utheta3];
    UFi_Cons = [C_uphi1;C_uphi2;C_uphi3];
    UPsi_Cons = [C_upsi1;C_upsi2;C_upsi3];

else 
    if Delta>0

        Q1 = DA*C_r1+3*((-DB+sqrt(Delta))/2);
        Q2 = DA*C_r1+3*((-DB-sqrt(Delta))/2);
        X = (-2*C_r1+(nthroot(Q1,3)+nthroot(Q2,3)))/6;
        Y = sqrt(3)*(nthroot(Q1,3)-nthroot(Q2,3))/6;
        alpha2 = sqrt(sqrt(X^2+Y^2)+X)/sqrt(2);
        beta2 = sqrt(sqrt(X^2+Y^2)+X)*(sqrt(X^2+Y^2)-X)/sqrt(2)/Y;
        gamma2 =sqrt((-C_r1-(nthroot(Q1,3)+nthroot(Q2,3)))/3);

        C_utheta4 = (alpha2^2-3*beta2^2)*C_theta1+C_theta2+C_theta3/(alpha2^2+beta2^2);
        C_utheta5 = (3*alpha2^2-beta2^2)*C_theta1+C_theta2-C_theta3/(alpha2^2+beta2^2);
        C_utheta6 = gamma2^3*C_theta1+gamma2*C_theta2+C_theta3/gamma2;
        C_uphi4 = C_phi1+(alpha2^2*C_utheta4-beta2^2*C_utheta4-2*beta2^2*C_utheta5)*C_phi2+C_utheta4*C_phi3;
        C_uphi5 = C_phi1+(2*alpha2^2*C_utheta4+alpha2^2*C_utheta5-beta2^2*C_utheta5)*C_phi2+C_utheta5*C_phi3;
        C_uphi6 = gamma2*C_phi1+gamma2^2*C_utheta6*C_phi2+C_utheta6*C_phi3;
        C_upsi4 = C_psi1+alpha2^2*C_utheta4*C_psi1-beta2^2*C_utheta5*C_psi1+alpha2^2*C_uphi4*C_psi2-beta2^2*C_uphi5*C_psi2;
        C_upsi5 = alpha2*beta2*(C_utheta4*C_psi1+C_utheta5*C_psi1+C_uphi4*C_psi2+C_uphi5*C_psi2);
        C_upsi6 = C_psi1+gamma2*C_utheta6*C_psi1+gamma2*C_uphi6*C_psi2;

        D0 = [1 1 gamma2;C_utheta4 C_utheta5 C_utheta6;C_uphi4 C_uphi5 C_uphi6];
        HD4 = [C_utheta4 C_utheta6;C_uphi4 C_uphi6];
        HD5 = [C_utheta5 C_utheta6;C_uphi5 C_uphi6];
        HD6 = [C_utheta4 C_utheta5;C_uphi4 C_uphi5];
        H4 = det(HD4)/det(D0);
        H5 = det(HD5)/det(D0);
        H6 = det(HD6)/det(D0);

        C1 = -Neutral_R*(alpha2*H4*cosh(alpha2*pi)*sin(beta2*pi)+beta2*H5*sinh(alpha2*pi)*cos(beta2*pi))/2/A55/alpha2/beta2/(cosh(alpha2*pi)^2*sin(beta2*pi)^2+sinh(alpha2*pi)^2*cos(beta2*pi)^2);
        C3 = Neutral_R*(alpha2*H4*sinh(alpha2*pi)*cos(beta2*pi)-beta2*H5*cosh(alpha2*pi)*sin(beta2*pi))/2/A55/alpha2/beta2/(cosh(alpha2*pi)^2*sin(beta2*pi)^2+sinh(alpha2*pi)^2*cos(beta2*pi)^2);
        C5 = -Neutral_R*H6/2/A55/sinh(gamma2*pi);

        U_Cons = [C1;C3;C5];
        Roots = [alpha2;beta2;gamma2];

        UTh_Cons = [C_utheta4;C_utheta5;C_utheta6];
        UFi_Cons = [C_uphi4;C_uphi5;C_uphi6];
        UPsi_Cons = [C_upsi4;C_upsi5;C_upsi6];

    end
end

end