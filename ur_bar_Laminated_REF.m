function [y] = ur_bar_Laminated_REF(x,rot,U_Cons,Roots,Delta)

C1 = U_Cons(1);
C3 = U_Cons(2);
C5 = U_Cons(3);

alpha = Roots(1);
beta = Roots(2);
gamma = Roots(3);

X = abs(x - rot);

if X > pi
    X = 2 * pi - X;

end

if Delta<0

    y = C1*cosh(alpha*X)+C3*cosh(beta*X)+C5*cosh(gamma*X);

else
    
    y = C1*cosh(alpha*X)*cos(beta*X)+C3*sinh(alpha*X)*sin(beta*X)+C5*cosh(gamma*X);

end

end