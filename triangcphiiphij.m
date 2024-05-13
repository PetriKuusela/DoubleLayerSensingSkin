function int=triangcphiiphij(g,c,ip,L, w)
%integrating with Gauss quadrature

phiatip = [1-ip(:,1)-ip(:,2),ip(:,1),ip(:,2)];
% Lasketaan c:n arvot integrointipisteissä
coupl = phiatip*c;

Jt = L*g;
dJt = abs(det(Jt));
bint = phiatip'*diag(coupl.*w)*phiatip;
int = dJt*bint;
