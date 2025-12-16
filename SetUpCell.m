function [Bsolar,Fsolar,nReal,nImag,N_EDphases,W_EDphases,Nlayers] = ...
    SetUpCell(materials,t,indActiveLayer,ActiveLayerVis,indLastOpticsLayer,lambda,AM15,opticalData,consts)
%SetUpPSC Assigns variables to workspace. 
% Inputs:
%   materials: material names
%   t: layer thicknesses (nm)
%   indActiveLayer: index of the active layer
%   ActiveLayerVis: interference visibility (/'coherence level') of the active layer
%   indLastOpticsLayer: index of the last material considered 
%   lambda: wavelengths (nm)
%   AM15: solar spectrum (nm,Wm^-2nm^-1)
%   opticalData: complex refractive indices of materials
%   consts: structure array of some natural constants

% Irradiation data
lambda2 = AM15{:,1};                        % Wavelength (nm)
B = AM15{:,2};                              % Radiance (Wm^-2nm^-1)

% Irradiance for defined wavelength range
Bsolar = spline(lambda2,B,lambda);
Fsolar = (lambda*1e-9)/(consts.h*consts.c).*Bsolar; % Photon flux (s^-1m^-2nm^-1)

% Useful variables
Nlayers = length(materials);

% Refractive indices and extinction coefficients of device materials
nReal = NaN(length(lambda),Nlayers+2);
nImag = NaN(length(lambda),Nlayers+2);

% Air before and after device
nReal(:,[1,end]) = ones(length(lambda),2);
nImag(:,[1,end]) = zeros(length(lambda),2);

% Interpolation
for i = 1:Nlayers
    if ~strcmp(materials{i},'Active')
        nReal(:,i+1) = interp1(rmmissing(opticalData.(strcat(materials{i},'_wl'))),...
            rmmissing(opticalData.(strcat(materials{i},'_n'))),lambda,'makima');
        nImag(:,i+1) = interp1(rmmissing(opticalData.(strcat(materials{i},'_wl'))),...
            rmmissing(opticalData.(strcat(materials{i},'_k'))),lambda,'makima');
    end
end

% Incoherent layers must be averaged, N_EDphases describes how many
% equidistant phases are used to calculate time average of specific layer
% (for coherent layers N_EDphases=1)
% Here, 1 micron is used as a threshold for coherent (<=1) and incoherent
% layers (>1 micron). In reality, coherency is continuous property and
% partial coherency should be applied.
N_EDphases = zeros(1,Nlayers);
for i = 1:Nlayers
    if t(i)<=1000
        N_EDphases(i) = 1;
    else
        N_EDphases(i) = 3;
    end
end

% Define weights for partial coherence in optical model
if ActiveLayerVis ~= 1
    N_EDphases(indActiveLayer) = 3;
end
M = prod(N_EDphases(1:indLastOpticsLayer));
Wlayers = cell(indLastOpticsLayer,1);
for i = 1:indLastOpticsLayer
    Wlayers{i} = 1/N_EDphases(i)*ones(N_EDphases(i),1);
end
W1PVK = (1-ActiveLayerVis)/N_EDphases(indActiveLayer)+ActiveLayerVis;
Wlayers{indActiveLayer} = [W1PVK;ones(N_EDphases(indActiveLayer)-1,1)*...
    (1-W1PVK)/(N_EDphases(indActiveLayer)-1)];
Ws = ones(prod(N_EDphases(1:indLastOpticsLayer)),indLastOpticsLayer);
count = 1;
for i = 1:indLastOpticsLayer
    if N_EDphases(i)~=1
        Ws(:,i) = repmat(repelem(Wlayers{i},count),[M/(count*length(Wlayers{i})),1]);
        count = count*N_EDphases(i);
    end
end
W_EDphases = prod(Ws,2);
end

