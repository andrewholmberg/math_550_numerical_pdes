% Constructing Gaussian Matrix

theta = linspace(0, 2*pi, 10);
theta(end) = [];  % Exclude the last element
xib = pi * cos(theta);  % Gaussian x-coordinate points
yib = pi * sin(theta);  % Gaussian y-coordinate points
N = 300;
xx = linspace(0, 2*pi, N);
yy = linspace(0, 2*pi, N);

dx = xx(2) - xx(1);

[X, Y] = meshgrid(xx, yy);  % Create a 2D grid
[Xint, Yint] = meshgrid(xx(2:end-1), yy(2:end-1));  % Interior points grid

% Define the function to spread values from the grid based on a delta

% Create the second derivative matrix (D2) using spdiags
e = (1 / dx^2) * ones(N, 1);
D2 = spdiags([e -2*e e], -1:1, N, N);

% Create identity matrix I_n
I_n = speye(N);

% Define the 2D Laplacian using Kronecker products
Lap = kron(I_n, D2) + kron(D2, I_n);

% LU decomposition (MATLAB's equivalent of scipy.sparse.linalg.splu)
dLap = decomposition(Lap, 'lu');



% Concatenation of vector for u and f
u = [reshape(cos(2*X) .* sin(2*Y), [], 1); zeros(length(xib), 1)];
f = [reshape(-8 * sin(2*Y) .* cos(2*X), [], 1); cos(2 * xib) .* sin(2 * yib)];



% Running gmres_helper function in MATLAB
xxx = gmres_helper(u, X, Y, xib, yib, @delta, Lap);

% gmres function to solve the system
[solution, flag] = gmres(@(v) gmres_helper(v, X, Y, xib, yib, @delta, Lap), f, [], 1e-6, 100);

disp('Solution:');
disp(solution);
disp('Flag:');
disp(flag);



% gmres_helper function
function output = gmres_helper(v, X, Y, xib, yib, delta, Lap)
    Sib = spreadQ(X, Y, xib, yib, delta);
    lq = length(xib);
    q = v(end-lq+1:end);
    u = v(1:end-lq);

    dx = X(1, 2) - X(1, 1);
    dy = Y(2, 1) - Y(1, 1);
    J = dx * dy * ones(lq, numel(X));
    tempJ = J * u;
    
    output = [Lap * u - Sib * q; tempJ];
end


function Sib = spreadQ(X, Y, xq, yq, delta)
    Nq = length(xq);
    Sib = zeros(Nq, numel(X));
    for k = 1:Nq
        Rk = sqrt((X - xq(k)).^2 + (Y - yq(k)).^2);
        Sib(k, :) = delta(Rk,:);  % Reshape into a column
    end
    Sib = Sib.';  % Transpose the result
end

% Define the function for interpolating Phi values based on delta
function Jphi = interpPhi(X, Y, xq, yq, Phi, delta)
    Nq = length(xq);
    Jphi = zeros(size(xq));
    dx = X(1, 2) - X(1, 1);
    dy = Y(2, 1) - Y(1, 1);
    
    for k = 1:Nq
        Rk = sqrt((X - xq(k)).^2 + (Y - yq(k)).^2);
        Jphi(k) = dx * dy * sum(sum(Phi .* delta(Rk)));
    end
end


% Define the Gaussian function delta
function delta_val = delta_a(r, a)
    delta_val = (1 / (2 * pi * a^2)) * exp(-0.5 * (r / a).^2);
end

function delta_val = delta(r)
    dx = 2*pi / 300;  % You need to ensure dx is defined properly
    delta_val = delta_a(r, 1.2 * dx);
end