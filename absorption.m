function [Gx,Glx,QxThermalization,Qthermalization,QxParasiticAbs,...
    QparasiticAbs,QxAbsorption,AbsLayers] = ...
    absorption(E2,x_pos,x_mat,n_real,k,EgPVK,EgPVKwl,...
    indPVKlayer,B_solar,lambda,Nlayers,consts,T)
%   E2(x,lambda): The square of electric field (E^2)
%   x_mat(1,x): Layer indices of x_pos
%   n_real(lambda,N_materials+2): real part of refractive indices for
%   materials and ambient media
%   k(lambda,N_materials+2): imaginary part of refractive indices 
%   (extinction coefficient)
%   Eg: Bandgap vector (0 if layer is not an active layer) (eV)
%   B_solar(lambda,1): Irradiance (W/(m^2nm)
%   lambda(lambda,1): Wavelengths (nm)

% Transpose E2
E2 = transpose(E2);

Qthermalization = zeros(length(lambda),length(x_pos));
QparasiticAbs = zeros(length(lambda),length(x_pos));
AbsLayers = zeros(length(lambda),Nlayers);
for matindex = 1:Nlayers
    % Position vector corresponding to current layer
    pos = x_mat==matindex;
    % Absorption coefficient in 1/m
    a=4*pi*k(:,matindex+1)./(lambda*1e-9);
    % Absorption
    Alx = repmat(a.*n_real(:,matindex+1),1,sum(pos)).*E2(:,pos);
    % Absorption in current layer for each wavelength
    AbsLayers(:,matindex) = trapz(x_pos(pos)*1e-9,Alx,2);
    % Absorbed energy
    AElx = Alx.*B_solar;
    if matindex == indPVKlayer
        % Wavelengths with enough energy to generate excitons
        Glambdas = lambda<=EgPVKwl;
        % Exciton generation rate
        Glx = AElx(Glambdas,:).*repmat((lambda(Glambdas)*1e-9)/...
            (consts.h*consts.c),1,sum(pos));
        Gx = zeros(1,sum(pos));
        for idxPos = 1:sum(pos)
            Gx(idxPos) = numInt_v1(lambda(1),EgPVKwl,transpose(lambda(Glambdas)),...
                transpose(Glx(Glambdas,idxPos)));
        end
        % Thermalization
        Qthermalization(Glambdas,pos) = Glx.*((consts.h*consts.c)./...
            (lambda(Glambdas)*1e-9)-(consts.q*EgPVK+3*consts.kB*T));
        QparasiticAbs(~Glambdas,pos) = AElx(~Glambdas,:);
    else
        QparasiticAbs(:,pos) = AElx;
    end
end

QxThermalization = trapz(lambda,Qthermalization,1);
QxParasiticAbs = trapz(lambda,QparasiticAbs,1);
% Q = Qthermalization + QparasiticAbs;
QxAbsorption = QxThermalization + QxParasiticAbs;
end

