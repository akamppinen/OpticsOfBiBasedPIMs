function I = Illuminance(B,effPhotopicVision)
% Illuminance calculates illuminance value for a given radiance
%Input:
% B: [wavelength (nm), radiance (W m^-2 nm^-1)] array
% effPhotopicVision: [wavelength (nm), V] array
%Output:
% I: illuminance (lux)

wl = 360:1:830; % nm
B = spline(B(:,1),B(:,2),wl);
effPhotopicVision = spline(effPhotopicVision(:,1),effPhotopicVision(:,2),wl);
I = 683*sum(B.*effPhotopicVision); % lux
end
