function [L_ex, ABD_ex]=ABD_curved_ex(Neutral_R,b,h,he,E_ex,nu_ex)

Q11 = E_ex*(1-nu_ex)/(1+nu_ex)/(1-2*nu_ex);
Q12 = E_ex*nu_ex/(1+nu_ex)/(1-2*nu_ex);

h_rho0 = (Neutral_R+h/2)*log((Neutral_R+h/2+he)/(Neutral_R+h/2));
h_rho1 = (Neutral_R+h/2)*(he - h_rho0);
h_rho2 = (Neutral_R+h/2)*he^2/2 - he*(Neutral_R+h/2)^2 + h_rho0*(Neutral_R+h/2)^2;

l1 = Q11*b*he;
l2 = Q11*b*he^2/2;
a11_e = Q12*b*h_rho0;
b11_e = Q12*b*h_rho1;
d11_e = Q12*b*h_rho2;

L_ex = [l1 l2];
ABD_ex = [a11_e b11_e d11_e];

end
