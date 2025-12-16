function [x_pos,x_mat,dx,t_cumsum] = xDiscretization(t,N,p)
%xDiscretization Creates one dimensional discretization
% Inputs
%   t (size: 1 x N of materials): Thicknesses of material layers
%   N (size: 1 x N of materials): Numbers of grid points in each layer
%   p (size: 1 x N of materials): Power factor determining grid point
%       density distributions in layers (distribution is denser towards the
%       the back side of material if p<0, uniform if p=0, and front dense 
%       if p>0.
% Outputs
%   x_pos (size: 1 x sum(N)): Grid points
%   x_mat (size: 1 x sum(N)): Corresponding layer numbers of grid points
%   dx (size: 1 x sum(N)): Corresponding slab thicknesses of grid points
%   t_cumsum (size: 1 x (N of materials + 1)): Cumulative sum of thicknesses
%   heatSourceBoundaries (size: 1 x (sum(N-1))): Slab boundary positions

% Cumulative sum of thicknesses
t_cumsum = cumsum([0,t]);

% Round N to integers
N = round(N);

% x_pos
x_pos = zeros(1,sum(N));
dx = zeros(1,sum(N));
% xElayerBoundaries = zeros(1,sum(N));
for material = 1:length(t)
    n = 1:N(material);
    mat = sum(N(1:material-1))+1:sum(N(1:material));
    dx(mat) = t(material)*n.^p(material)/sum(n.^p(material));
    x_pos(mat) = cumsum([t_cumsum(material),dx(mat(1:end-1))])+dx(mat)/2;
%     xElayerBoundaries(mat) = ...
%         t_cumsum(material)+cumsum(dx(mat(1:end)));
end

% x_mat
x_mat = sum(repmat(x_pos,length(t_cumsum),1)>repmat(t_cumsum',1,length(x_pos)),1);
end