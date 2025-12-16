function [E,M] = SystemMatrix(n,k,d,lambda,N_EDphases)
%SystemMatrix constructs System matrix (set of equations of light
%propagation in given layers) and solves electric field amplitudes at the
%layer interfaces.
% n: real part of refractive index
% k: imaginary part of refractive index
% d: layer thicknesses (nm)
% lambda: wavelengths (nm)

I = size(n,2)-1;
J = size(n,1);
M = prod(N_EDphases);

E = zeros(J,4*I,M);
t = zeros(J,I);
r = zeros(J,I);
delta = zeros(J,I);

% Phase shifts for averaging incoherent layers
phi = zeros(M,length(N_EDphases));
count = 1;
for i = 1:length(N_EDphases)
    if N_EDphases(i)~=1
        EDs = linspace(0,(1-1/N_EDphases(i))*2*pi,N_EDphases(i))';
        phi(:,i) = repmat(repelem(EDs,count),[M/(count*length(EDs)),1]);
        count = count*N_EDphases(i);
    end
end

for i = 1:I
    % Reflection and transmission coefficients for all interfaces from Fresnel 
    % equations
    [t(:,i),r(:,i),~,~,~] = Fresnel_complex_Vector(n(:,i),k(:,i),n(:,i+1),k(:,i+1),0);
    
    % Calculate phase changes for coherent light propagation (last "layer"
    % is considered infinite, and it is typically air)
    if i<I
        delta(:,i) = 2*pi*(n(:,i+1)-k(:,i+1)*1i)*d(i)./lambda;
    end
end

y = zeros(4*I,1); y(1) = 1;
% For all wavelengths
for j = 1:J
    for m = 1:M
    % Construct System matrix A
    A = eye(4*I);
    for i = 1:I 
        base = (i-1)*4+1;
        if i == 1
            A(base+1,base+2) = -1/t(j,i);
            A(base+1,base+3) = -r(j,i)/t(j,i);
            A(base+2,base+5) = -exp(-(delta(j,i)+phi(m,i))*1i);
            A(base+3,base) = -t(j,i);
            A(base+3,base+2) = r(j,i);
        elseif i == I
            A(base,base-1) = -exp(-(delta(j,i-1)+phi(m,i-1))*1i);
            A(base+1,base+2) = -1/t(j,i);
            A(base+1,base+3) = -r(j,i)/t(j,i);
            A(base+3,base) = -t(j,i);
            A(base+3,base+2) = r(j,i);
        else
            A(base,base-1) = -exp(-(delta(j,i-1)+phi(m,i-1))*1i);
            A(base+1,base+2) = -1/t(j,i);
            A(base+1,base+3) = -r(j,i)/t(j,i);
            A(base+2,base+5) = -exp(-(delta(j,i)+phi(m,i))*1i);
            A(base+3,base) = -t(j,i);
            A(base+3,base+2) = r(j,i);
        end
    end
    % Solve E field amplitudes
    E(j,:,m) = A\y;
    end
end
assignin('base','CoeffMatrixA',A);
end