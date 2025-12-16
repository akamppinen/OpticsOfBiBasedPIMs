function [t_s,r_s,t_p,r_p,theta_t] = Fresnel_complex_Vector(n1,k1,n2,k2,theta_i)

% n1,n2: Real part of refractive indices of materials
% k1,k2: Imaginary part of refractive indices of materials
% theta_i: Incidence angle (degrees)

% Complex refractive indices
n1 = n1-k1*1i;
n2 = n2-k2*1i;

% Unit conversion
theta_i = 2*pi*theta_i./360;
% Angle of transmitted light based on Snell's law
theta_t = asin((n1./n2).*sin(theta_i));

% s polarization
t_s = 2*n1.*cos(theta_i)./(n1.*cos(theta_i)+n2.*cos(theta_t));
r_s = (n1.*cos(theta_i)-n2.*cos(theta_t))./(n1.*cos(theta_i)+n2.*cos(theta_t));
% p polarization
t_p = 2*n1.*cos(theta_i)./(n1.*cos(theta_t)+n2.*cos(theta_i));
r_p = (n1.*cos(theta_t)-n2.*cos(theta_i))./(n1.*cos(theta_t)+n2.*cos(theta_i));

% Unit conversion
theta_t = 360*theta_t./(2*pi);
end