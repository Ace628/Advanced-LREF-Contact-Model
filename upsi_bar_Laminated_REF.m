function [y] = upsi_bar_Laminated_REF(x,rot,U_Cons,Roots,Delta,UPsi_Cons)

C1 = U_Cons(1);
C3 = U_Cons(2);
C5 = U_Cons(3);

alpha = Roots(1);
beta = Roots(2);
gamma = Roots(3);

C_upsi1 = UPsi_Cons(1);
C_upsi2 = UPsi_Cons(2);
C_upsi3 = UPsi_Cons(3);

X = abs(x - rot);

if X > pi
    X = 2 * pi - X;
end

if Delta<0

    y = C1*C_upsi1*cosh(alpha*X)+C3*C_upsi2*cosh(beta*X)+C5*C_upsi3*cosh(gamma*X);

else

    y = C1*(C_upsi1*cosh(alpha*X)*cos(beta*X)-C_upsi2*sinh(alpha*X)*sin(beta*X))+...
        C3*(C_upsi1*sinh(alpha*X)*sin(beta*X)+C_upsi2*cosh(alpha*X)*cos(beta*X))+...
        C5*C_upsi3*cosh(gamma*X);

end

end