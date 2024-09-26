
[x,y,u,v] = solve_linear_poisson(100)
u
function [x_vals, y_vals, u_vals, v_vals] = solve_linear_poisson(N)
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
        f(n) = -8 * sin(2 * xx(j)) * cos(2 * xx(i));
        
        % Store x and y coordinates
        x_list(end+1) = (i - 1) * dx;
        y_list(end+1) = (j - 1) * dy;
        
        % Apply boundary conditions
        if i == 1 || i == N || j == 1 || j == N
            % Dirichlet boundary condition
            data(end+1) = 1;
            row_ind(end+1) = n;
            col_ind(end+1) = n;
            f(n) = u_exact(xx(i), xx(j));  % Set boundary value to exact solution
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

    % Solve for the numerical solution (flattened vector)
    u_n = A \ f;

    % Reshape the solution back into a 2D grid
    x_vals = reshape(x_list, N, N);
    y_vals = reshape(y_list, N, N);
    u_vals = reshape(u_n, N, N);
    
    % Compute the analytical solution on the same grid
    v_vals = u_exact(x_vals, y_vals);

    % Plot the numerical solution
    figure;
    surf(x_vals, y_vals, u_vals);
    xlabel('x');
    ylabel('y');
    zlabel('u(x, y)');
    title('Numerical Solution to the Poisson Equation');

    % Plot the analytical solution
    figure;
    surf(x_vals, y_vals, v_vals);
    xlabel('x');
    ylabel('y');
    zlabel('Exact u(x, y)');
    title('Analytical Solution to the Poisson Equation');
end

