% Constructing Gaussian Matrix
clear all;

N = 100;
m = 6;
theta = linspace(0, 2*pi, m);
theta(end) = [];  % Exclude the last element
xib = pi * cos(theta);  % Gaussian x-coordinate points
yib = pi * sin(theta);  % Gaussian y-coordinate points
xx = linspace(0, 2*pi, N);
yy = linspace(0, 2*pi, N);

dx = xx(2) - xx(1);

[X, Y] = meshgrid(xx, yy);  % Create a 2D grid
xvals = X(:);
yvals = Y(:);
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
p1 = -8 * sin(2.*Y) .* cos(2.*X);
p1 = p1(:);
p2 = cos(2 * xib) .* sin(2 * yib);
p2 = p2(:);
f = [p1; p2];
% Initialize the right-hand side (RHS) vector f (source term)

[A, f] = solve_linear_poisson(N);
size(A)
size(f)
u_n = A \ f;
f = [f;p2];
u = [u; zeros(length(xib), 1)];

% gmre
% ?""?s function to solve the system
[solution, flag] = gmres(@(v) gmres_helper(v, X, Y, xib, yib, @delta, A), f, [], 1e-6, 100);



u = solution(1:end-length(xib));
k = solution(end-length(xib)+1:end);


shp = sqrt(length(u));
U = reshape(u,shp,shp);
X = reshape(X,shp,shp);
Y = reshape(Y,shp,shp);
V = cos(2*X).*sin(2*Y);

[x_new, y_new] = meshgrid(0:0.1:2*pi, 0:0.1:2*pi); % New smaller grid size
% Interpolate data onto the new grid
z_new = interp2(X, Y, U, x_new, y_new, 'linear'); % 'linear' interpolation
% Display the loaded variables
% Create a 3D surface plot
figure; % Open a new figure window
surf(x_new, y_new, z_new); % Create a 3D surface plot
% Set axis labels
xlabel('X');
ylabel('Y');
zlabel('U');
% Add a title to the plot
title("Numerical Solution to Poisson's Equation Where u(x,y) = sin(x)cos(y)");
% Add a colorbar for reference
colorbar;
% Set the view angle for better visualization
view(45, 30); % Adjust azimuth and elevation for a better view
[x_new, y_new] = meshgrid(0:0.1:2*pi, 0:0.1:2*pi); % New smaller grid size
% Interpolate data onto the new grid

v_new = sin(2*y_new).*cos(2*x_new);
%v_new = interp2(X, Y, V, x_new, y_new, 'linear'); % 'linear' interpolation
% Display the loaded variables
% Create a 3D surface plot
figure; % Open a new figure window
surf(x_new, y_new, v_new); % Create a 3D surface plot
% Set axis labels
xlabel('X');
ylabel('Y');
zlabel('U');
% Add a title to the plot
title("Analytical Solution to Poisson's Equation Where u(x,y) = sin(x)cos(y)");
% Add a colorbar for reference
colorbar;
% Set the view angle for better visualization
view(45, 30); % Adjust azimuth and elevation for a better view


%}


% gmres_helper function
function output = gmres_helper(v,X,Y,xib,yib,delta,Lap)

    Sib = spreadQ(X, Y, xib, yib, delta);
    lq = length(xib);
    q = v(end-lq+1:end);
    u = v(1:end-lq);
    dx = X(1, 2) - X(1, 1);
    dy = Y(2, 1) - Y(1, 1);
    %J = dx * dy * ones(lq, numel(X));
    J = interpPhi(X,Y,xib,yib,delta);
    shape(J)
    output = [Lap * u - Sib * q; J * u];
end


function Sib = spreadQ(X, Y, xq, yq, delta)
    Nq = length(xq);
    Sib = zeros(Nq, numel(X));
    size(Sib)
    for k = 1:Nq
        Rk = sqrt((X - xq(k)).^2 + (Y - yq(k)).^2);
        x = reshape(delta(Rk),1,[]);
        Sib(k,:) = x.' ;  % Reshape into a column
    end
    Sib = Sib.';  % Transpose the result
end

% Define the function for interpolating Phi values based on delta
function Jphi = interpPhi(X, Y, xq, yq, delta)
    Nq = length(xq);
    Jphi = zeros(size(xq));
    dx = X(1, 2) - X(1, 1);
    dy = Y(2, 1) - Y(1, 1);
    
    for k = 1:Nq
        Rk = sqrt((X - xq(k)).^2 + (Y - yq(k)).^2);
        %Jphi(k) = dx * dy * sum(sum(Phi .* delta(Rk)));
        temp = reshape(delta(Rk),[],1);
        Jphi(k) = dx * dy * temp;
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





function [A,f] = solve_linear_poisson(N)
    % Define the grid
    xx = linspace(0, 2*pi, N);
    dx = xx(2) - xx(1);
    dy = dx;  % Assuming dx = dy

    % Initialize storage for sparse matrix construction
    data = [];
    row_ind = [];
    col_ind = [];
    x_list = [];
    y_list = [];
    
    % Function for exact solution (Dirichlet boundary conditions)
    u_exact = @(x, y) sin(2 * y) .* cos(2 * x);
    
    % Initialize right-hand side vector (source term)
    f = zeros(N * N, 1);  
    
    % Build the matrix and RHS vector
    for n = 1:N*N
        i = mod(n-1, N) + 1;  % Row index
        j = floor((n-1) / N) + 1;  % Column index
        
        % Compute f at each grid point
        f(n) = -8 * sin(2 * xx(i)) * cos(2 * xx(j));
        
        % Store x and y coordinates
        x_list(end+1) = (j - 1) * dx;
        y_list(end+1) = (i - 1) * dy;
        
        % Apply boundary conditions
        if i == 1 || i == N || j == 1 || j == N
            % Dirichlet boundary condition
            data(end+1) = 1;
            row_ind(end+1) = n;
            col_ind(end+1) = n;
            f(n) = u_exact(xx(j), xx(i));  % Set boundary value to exact solution
        else
            % Interior points: 5-point stencil for Laplacian
            neighbors = [-1, 0, 1, 0, 0];
            offsets = [-1, 0, 1, -N, N];  % Left, center, right, down, up

            for k = 1:length(neighbors)
                coef =  (1/dx^2);
                if k == 2  % center point
                    coef = -4 / dx^2;
                end
                % Add values to the sparse matrix
                row_ind(end+1) = n;
                col_ind(end+1) = n + offsets(k);
                data(end+1) = coef;
            end
        end
    end

    % Construct the sparse matrix
    A = sparse(row_ind, col_ind, data, N*N, N*N);
end