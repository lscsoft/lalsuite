function f = fofMa(M,a)

G = 6.6742e-11;        % m^3/kg/s^2
c = 299792458.;        % m/s
Msun = 1.9889e30;      % kg

f = 1/(2*pi) .* c^3./(G.*M*Msun) .* (1-0.63.*(1-a).^(3/10));