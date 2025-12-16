function [E2ave,E,R,A,T,Pave] = CalcSolarCellOptics_v4(n_real,k,t,x_pos,x_mat,lambda,N_EDphases,W_EDphases)
%CalcSolarCellOptics Solves the square of electric field as a function of
%position in a solar cell defined as input variables. In addition,
%reflection from the cell, absorption in different layers and transmission
%through the cell are calculated. A form of the transfer matrix method is
%applied here.
% Inputs
%   n_real:
%   k:
%   t:
%   x_pos:
%   x_mat:
%   lambda:
%   N_EDphases:
% Outputs
%   E2ave:
%   R:
%   A:
%   T:
%   Pave:

% Version 4.1
% Add weights for averaging equidistant phases

% Check variable dimensions
if size(lambda,2)~=1
    error('lambda must be a column vector.')
elseif size(t,1) ~= 1
    error('t must be a row vector.')
elseif size(x_pos,1) ~= 1
    error('x_pos must be a row vector.')
elseif size(x_mat) ~= size(x_pos)
    error('x_mat and x_pos must be of same size.')
elseif ~all(size(n_real) == [length(lambda),length(t)+2])
    error('n and k must be of size (N of wavelengths) x (N of layers+2).')
elseif size(n_real) ~= size(k)
    error('n and k must be of same size.')
elseif length(W_EDphases) ~= prod(N_EDphases)
    error('Give weights for all equidistant phases.')
elseif abs(sum(W_EDphases)-1)>1e-5
    error('Weights must sum up to one.')
end

% Cumulative sum of layers
t_cumsum = cumsum([0,t]);

% Solve E at material interfaces
[E_atInterfaces,M] = SystemMatrix(n_real,k,t,lambda,N_EDphases);

% Weights of different phases
Ws = permute(repmat(W_EDphases,[1,length(x_pos),length(lambda)]),[2,3,1]);

% E field in the layers
E = zeros(length(x_pos),length(lambda),M); % electric field
H = zeros(length(x_pos),length(lambda),M); % magnetic field
P = zeros(length(x_pos),length(lambda),M); % Poynting vector
for l = 1:length(lambda)
    for m = 1:M
        for material=1:length(t)%size(layers,2)
            base_column = (material-1)*4+1;
            E_forward = E_atInterfaces(l,base_column+3,m);
            E_backward = E_atInterfaces(l,base_column+5,m);
            
            xi = 2*pi*(n_real(l,material+1)-1i*k(l,material+1))/lambda(l);
            dj = t(material);
            %indices of points which are in the material layer
            x_indices = find(x_mat == material);
            %distance from interface with previous layer
            x = x_pos(x_indices)-t_cumsum(material); 
            
            E(x_indices,l,m) = (E_forward*exp(-1i*(xi*x))+...
                E_backward*exp(-1i*(xi*(dj-x))))/sqrt(n_real(l,1)); % with normalization!
            H(x_indices,l,m) = (n_real(l,material+1)-1i*k(l,material+1))./...
                (n_real(l,1)-1i*k(l,1)).*...
                (E_forward*exp(-1i*(xi*x))-E_backward*exp(-1i*(xi*(dj-x))));
            P(x_indices,l,m) = real(E(x_indices,l,m).*conj(H(x_indices,l,m)));
        end
    end
end
% Average Poynting vector 
Pave = sum(P.*Ws,3);
% Field intensities in x_positions for wavelengths
E2 = abs(E).^2;
% Average over the runs of equidistant phase shifts in incoherent layers
E2ave = sum(E2.*Ws,3);

% Poynting vector at the interfaces and absorption in layers 
E_tot = zeros(size(n_real,1),size(n_real,2)-1,M);
H_tot = zeros(size(n_real,1),size(n_real,2)-1,M);
P_atInterfaces = zeros(size(n_real,1),size(n_real,2)-1);
for m = 1:M
for i = 1:(size(n_real,2)-1)
    base_column = (i-1)*4+1;
    E_tot(:,i,m) = E_atInterfaces(:,base_column,m)+E_atInterfaces(:,base_column+1,m);
    H_tot(:,i,m) = (n_real(:,i)-k(:,i)*1i)./(n_real(:,1)-k(:,1)*1i).*...
        (E_atInterfaces(:,base_column,m)-E_atInterfaces(:,base_column+1,m));
    
    P_atInterfaces(:,i) = P_atInterfaces(:,i) + ...
        W_EDphases(m)*real(E_tot(:,i,m).*conj(H_tot(:,i,m)));
end
end
% P_atInterfaces = P_atInterfaces/M;

R = ones(length(lambda),1)-P_atInterfaces(:,1);
A = P_atInterfaces(:,1:(end-1))-P_atInterfaces(:,2:end);
T = P_atInterfaces(:,end);
end

