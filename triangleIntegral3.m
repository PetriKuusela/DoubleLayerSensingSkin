function int=triangleIntegral3(g, f1, f2, f3, ip, L, w)
%integrating f1*f2*f3 with Gauss quadrature.


phiatip = [1-ip(:,1)-ip(:,2),ip(:,1),ip(:,2)];
% Lasketaan funktioden arvot integrointipisteissä
fp1 = phiatip*f1;
fp2 = phiatip*f2;
fp3 = phiatip*f3;

Jt = L*g;
dJt = abs(det(Jt));
bint = sum(fp1.*fp2.*fp3.*w);
int = dJt*bint;
