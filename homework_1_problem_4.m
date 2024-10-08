% Constructing Gaussian Matrix
clear all;


u =@(x,y) cos(2.*x).*sin(2.*y);
f =@(x,y) -8*cos(2.*x).*sin(2.*y);
f_ib =@(x,y) cos(2.*x).*sin(2.*y);

% params
N = 500;
m=5;
tol = 10^-2;%(-1.5); % tolerance for gmres

% make grid domain
xdom = linspace(0,2*pi,N);
ydom = xdom;
h = xdom(2)-xdom(1);
xdom_inner = xdom(2:end-1);
ydom_inner = xdom_inner;

[Xdom,Ydom] = meshgrid(xdom(2:end-1),ydom(2:end-1));
[Xdom_full,Ydom_full] = meshgrid(xdom,ydom);
fmat = f(Xdom,Ydom);
u_exact = u(Xdom_full,Ydom_full);

% Add BCs (Dirichlet for grid)
fmat(1,:) = fmat(1,:) - (1/h^2)*u_exact(1,2:end-1);
fmat(end,:) = fmat(end,:) - (1/h^2)*u_exact(end,2:end-1);
fmat(:,1) = fmat(:,1) - (1/h^2)*u_exact(2:end-1,1);
fmat(:,end) = fmat(:,end) - (1/h^2)*u_exact(2:end-1,end);

% make IB domain
% N_ib = 100;
theta_all = 0:h:2*pi;%linspace(0,2*pi,N_ib); % careful, makes 1 pt overlap
theta = theta_all(1:end-1); % avoids overlap of theta=0,2pi
N_ib = length(theta);
xib = pi + cos(theta);
yib = pi + sin(theta);
% [Xdom_ib,Ydom_ib] = meshgrid(xib(1:end),yib(1:end));
% [Xdom_ib_full,Ydom_ib_full] = meshgrid(xib,yib); % not needed
% h_ib = ***
% fmat_ib = f_ib(Xdom_ib,Ydom_ib); % matrix of immersed boundary conditions for fvector
fmat_ib = f_ib(xib,yib); % vector of immersed boundary conditions for fvector

% delta function approx wt gaussian
delta_a = @(r,a) (1/(2*pi*a^2))*exp(-0.5*(r/a).^2); 
delta = @(r) delta_a(r,1.2*h);

% Reshape f matrixes to vector
% fvect = [reshape(fmat,(N-2)^2,1); reshape(fmat_ib,(N_ib),1)];
fvect = [reshape(fmat,[],1)*0; fmat_ib'];

% FDM matrix
evect = ones(N-2,1);
D2 = (1/h^2)*spdiags([evect -2*evect evect], -1:1, N-2, N-2);
I_n = speye(N-2);
D_lap = kron(I_n, D2) + kron(D2, I_n);

% Solve for unknowns 
% afun =@(x) afun_N(x,N,N_ib,X,Y,x_ib,y_ib,delta,D_lap); %recast of funct only of uvect
% u_appvect = gmres(@afun, fvect, 1000, tol, 1); % solver
Xdom_vect = reshape(Xdom,(N-2)^2,1);
Ydom_vect = reshape(Ydom,(N-2)^2,1);

% uu = ones((N-2)^2+N_ib,1);
% test_fun = afun_N(uu,N,N_ib,Xdom_vect,Ydom_vect,xib,yib,delta,D_lap,h);
indices = return_points_indices(h,min(Xdom_vect),max(Xdom_vect),min(Ydom_vect),max(Ydom_vect),xib,yib,.0001);

helper = @(x)gmres_helper(x,N,N_ib,Xdom_vect,Ydom_vect,xib,yib,delta,D_lap,h,indices);
tic
[solution, flag] = gmres(helper, fvect, [], tol, 1000); % solver
toc

u_solution = solution(1:end-length(xib));
k = solution(end-length(xib)+1:end);

shp = sqrt(length(u_solution));
U = reshape(u_solution,shp,shp);
U = [u_exact(1,2:end-1)*0; U; u_exact(end,2:end-1)*0];
U = [u_exact(:,1)*0 U u_exact(:,end)*0];
X = reshape(Xdom_full,shp+2,shp+2);
Y = reshape(Ydom_full,shp+2,shp+2);
V = cos(2*X).*sin(2*Y);
uib = u(xib,yib);
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
plot3(xib,yib,uib,'ko','MarkerSize',6,'MarkerFaceColor', 'w')
hold on
%colormap('gray')
surf(x_new, y_new, z_new); % Create a 3D surface plot

grid on;% Set axis labels
xlabel('X');
ylabel('Y');
zlabel('U');
% Add a title to the plot
title("Analytical Solution to Poisson's Equation Where u(x,y) = sin(x)cos(y)");
% Add a colorbar for reference
colorbar;
% Set the view angle for better visualization
view(45, 30); % Adjust azimuth and elevation for a better view



%% Functions
% mat to solve system ***
function [Sq] = spreadQ(X,Y,xib,yib,N_ib,q,delta,indices)
    Sq = 0*X;
    % Nq = length(q);
    for k = 1:N_ib
        Rk = sqrt((X(indices(k))-xib(k)).^2 + (Y(indices(k))-yib(k)).^2);
        Sq(indices(k)) = Sq(indices(k)) + q(k)*delta(Rk);
    end
end

function [Jphi] = interpPhi(X,Y,xib,yib,N_ib,Phi,delta,h,indices)
    Jphi = 0*xib;
    dx = h;
    dy = h;

    for k = 1:N_ib
        Rk = sqrt((X(indices(k))-xib(k)).^2 + (Y(indices(k))-yib(k)).^2);
        Jphi(k) = dx*dy*sum((Phi(indices(k)).*delta(Rk)));
        % Jphi(k) = dx*dy*sum(sum(Phi.*delta(Rk)));
    end
end

function A = gmres_helper(x,N,N_ib,X,Y,xib,yib,delta,D_lap,h,indices)
    vec_len = (N-2)^2 + N_ib;
    u_len = (N-2)^2;

    % apply functions to u
    Ax(1:u_len) = D_lap*x(1:u_len);
    Ax(u_len+1:vec_len) = interpPhi(X,Y,xib,yib,N_ib,x(1:u_len),delta,h,indices); %J(x(1:N-2)),  J(u)

    % apply functions to q
    % test = spreadQ(X,Y,xib,yib,N_ib,x(u_len+1:vec_len),delta) ;
    Ax(1:u_len) =  Ax(1:u_len) - (spreadQ(X,Y,xib,yib,N_ib,x(u_len+1:vec_len),delta,indices))' ; %Ax(1:N-2) - S(x(N-1:vec_len)); -S(q)
    
    A = Ax';
end

% afun =@(x) afun_N(x,N,N_ib,X,Y,x_ib,y_ib,delta,D_lap); %recast of funct only of uvect

% function out = afun(x) 
%     out = afun_N(x,N,N_ib,X,Y,x_ib,y_ib,delta,D_lap); %recast of funct only of uvect
% end
function m =  return_points_indices(h,xmin,xmax,ymin,ymax,xib,yib,tol)
    reach = abs(norminv(tol,0,1)) * h;
    nxpoints = floor((xmax-xmin)/h)+1;
    nypoints = floor((ymax-ymin)/h)+1;
    m = containers.Map([1], {[1 8;]});
    for i = 1:length(xib)
        xmind = floor((xib(i) - reach-xmin)/h);
        xmaxd = ceil((xib(i) + reach-xmin)/h);
        ymind = floor((yib(i) - reach-ymin)/h);
        ymaxd =  ceil((yib(i) + reach-ymin)/h);
        indices = zeros((xmaxd-xmind) * (ymaxd-ymind));
        indices = [];
        for x = xmind : xmaxd+1
            indices = [indices; (ymind + x * nypoints):1:(ymaxd + x * nypoints)];
        end
        
        %[xx,yy] = meshgrid(xmind:1:xmaxd,ymind:1:ymaxd);
        %xx=reshape(xx,[],1);
        %yy=reshape(yy,[],1);
        indices = reshape(indices,[],1);
        indices = indices(indices >0 & indices <= nxpoints * nypoints);
        m(i) = unique(indices);

    end


    
end 