function [y] = ut_bar_Laminated_REF(x,rot,U_Cons,Roots,Delta,UTh_Cons)

C1 = U_Cons(1);
C3 = U_Cons(2);
C5 = U_Cons(3);

alpha = Roots(1);
beta = Roots(2);
gamma = Roots(3);

C_utheta1 = UTh_Cons(1);
C_utheta2 = UTh_Cons(2);
C_utheta3 = UTh_Cons(3);

X = abs(x - rot);

if X > pi
    X = 2 * pi - X;
end

if Delta<0

    y = C1*C_utheta1*sinh(alpha*X)+C3*C_utheta2*sinh(beta*X)+C5*C_utheta3*sinh(gamma*X);

else

    y = C1*(alpha*C_utheta1*sinh(alpha*X)*cos(beta*X)-beta*C_utheta2*cosh(alpha*X)*sin(beta*X))+...
        C3*(alpha*C_utheta1*cosh(alpha*X)*sin(beta*X)+beta*C_utheta2*sinh(alpha*X)*cos(beta*X))+...
        C5*C_utheta3*sinh(gamma*X);

end

if rot > 0
    if x>=-pi+rot && x<= rot
        y = -y;
    end
elseif rot < 0
    if x>=-pi && x<=rot
        y = -y;
    elseif x >= pi+rot && x<=pi
        y = -y;
    end
elseif rot == 0
    if x>=-pi && x<=0
        y = -y;
    end
end

end