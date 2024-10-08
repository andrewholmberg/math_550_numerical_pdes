%% MATH 550 HW 1 P2a
clear all;
close all;

fontsize = 18;
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(0, 'defaultTextInterpreter', 'latex');
set(0, 'DefaultFigureRenderer', 'painters');
set(groot,'defaultAxesFontSize',fontsize);

%% Solve Problem
% parameter functions
u =@(x,y) cos(2.*x).*sin(2.*y);
f =@(x,y) -8*cos(2.*x).*sin(2.*y);
f_ib =@(x,y) cos(2.*x).*sin(2.*y);

% params
N = 10;
tol = 10^-6;%(-1.5); % tolerance for gmres

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
x_ib_dom = pi + cos(theta);
y_ib_dom = pi + sin(theta);
% [Xdom_ib,Ydom_ib] = meshgrid(x_ib_dom(1:end),y_ib_dom(1:end));
% [Xdom_ib_full,Ydom_ib_full] = meshgrid(x_ib_dom,y_ib_dom); % not needed
% h_ib = ***
% fmat_ib = f_ib(Xdom_ib,Ydom_ib); % matrix of immersed boundary conditions for fvector
fmat_ib = f_ib(x_ib_dom,y_ib_dom); % vector of immersed boundary conditions for fvector

% delta function approx wt gaussian
delta_a = @(r,a) (1/(2*pi*a^2))*exp(-0.5*(r/a).^2); 
delta = @(r) delta_a(r,1.2*h);

% Reshape f matrixes to vector
% fvect = [reshape(fmat,(N-2)^2,1); reshape(fmat_ib,(N_ib),1)];
fvect = [reshape(fmat,(N-2)^2,1); fmat_ib'];

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
% test_fun = afun_N(uu,N,N_ib,Xdom_vect,Ydom_vect,x_ib_dom,y_ib_dom,delta,D_lap,h);

afun = @(x)afun_N(x,N,N_ib,Xdom_vect,Ydom_vect,x_ib_dom,y_ib_dom,delta,D_lap,h);
u_appvect = gmres(afun, fvect, [], tol, 1000); % solver
tester = spreadQ(Xdom_vect,Ydom_vect,x_ib_dom,y_ib_dom,N_ib,ones(size(x_ib_dom)),delta)
% u_appvect = gmres(@(x)afun_N(x,N,N_ib,Xdom,Ydom,x_ib_dom,y_ib_dom,delta,D_lap), fvect, [], tol, 1000); % solver
% u_appvect = A\fvect;
u_approx = reshape(u_appvect(1:(N-2)^2),N-2,N-2);
% q_approx = reshape(u_appvect((N-2)^2+1:end), N_ib,1);
q_approx = u_appvect((N-2)^2+1:end);

% Interp interior BC ***
Phi_ib = interpPhi(Xdom_vect,Ydom_vect,x_ib_dom,y_ib_dom,N_ib,u_appvect(1:(N-2)^2),delta,h);
% Phi_ib = reshape(Phi_ib_mat,N-2,N-2);
% interpPhi(Xdom_vect,Ydom_vect,xib_dom,yib_dom,N_ib,Phi,delta,h)

% Add BCs to u (phi)
u_approx = [u_exact(1,2:end-1); u_approx; u_exact(end,2:end-1)];
u_approx = [u_exact(:,1) u_approx u_exact(:,end)];

% error
norm_type =2;
rel_error_2 = norm(u_exact-u_approx,norm_type)/norm(u_exact,norm_type)

%% Plotting
figure()
colormap('cool')
mesh(Xdom_full,Ydom_full,u_approx,'FaceAlpha', 0.7,'EdgeAlpha',0.8)
hold on
% plot3(x_ib_dom,y_ib_dom,q_approx,'or')
% plot3(x_ib_dom,y_ib_dom,fmat_ib,'.b', 'MarkerSize',15)
plot3(x_ib_dom,y_ib_dom,Phi_ib,'.-b', 'MarkerSize',20)
xlabel('$x$')
ylabel('$y$')
zlabel('$\Phi_{approx}(x,y)$')
title('Approximate Solution')
% zlim([-1,1])
legend('$\Phi$', 'Imm. Bound.')

figure()
mesh(Xdom_full,Ydom_full,u_exact)
hold on
% plot3(x_ib_dom,y_ib_dom,fmat_ib,'ob')
colormap('cool')
xlabel('$x$')
ylabel('$y$')
zlabel('$u_{exact}(x,y)$')
title('Exact Solution')

% figure()
% hold on
% plot3(x_ib_dom,y_ib_dom,q_approx,'ob')
% plot3(x_ib_dom,y_ib_dom,f_ib(x_ib_dom,y_ib_dom),'om')

%% Functions
% mat to solve system ***
function [Sq] = spreadQ(X,Y,xib,yib,N_ib,q,delta)
    Sq = 0*X;
    % Nq = length(q);
    for k = 1:N_ib
        Rk = sqrt((X-xib(k)).^2 + (Y-yib(k)).^2);
        Sq = Sq + q(k)*delta(Rk);
    end
end

function [Jphi] = interpPhi(X,Y,xib,yib,N_ib,Phi,delta,h)
    Jphi = 0*xib;
    dx = h;
    dy = h;
    for k = 1:N_ib
        Rk = sqrt((X-xib(k)).^2 + (Y-yib(k)).^2);
        Jphi(k) = dx*dy*sum((Phi.*delta(Rk)));
        % Jphi(k) = dx*dy*sum(sum(Phi.*delta(Rk)));
    end
end

function A = afun_N(x,N,N_ib,X,Y,xib,yib,delta,D_lap,h)
    vec_len = (N-2)^2 + N_ib;
    u_len = (N-2)^2;

    % apply functions to u
    Ax(1:u_len) = D_lap*x(1:u_len);
    Ax(u_len+1:vec_len) = interpPhi(X,Y,xib,yib,N_ib,x(1:u_len),delta,h); %J(x(1:N-2)),  J(u)

    % apply functions to q
    % test = spreadQ(X,Y,xib,yib,N_ib,x(u_len+1:vec_len),delta) ;
    Ax(1:u_len) =  Ax(1:u_len) - (spreadQ(X,Y,xib,yib,N_ib,x(u_len+1:vec_len),delta))' ; %Ax(1:N-2) - S(x(N-1:vec_len)); -S(q)
    
    A = Ax';
end

% afun =@(x) afun_N(x,N,N_ib,X,Y,x_ib,y_ib,delta,D_lap); %recast of funct only of uvect

% function out = afun(x) 
%     out = afun_N(x,N,N_ib,X,Y,x_ib,y_ib,delta,D_lap); %recast of funct only of uvect
% end


%% Error Loop
% clear;
% 
% u =@(x,y) sin(x).*cos(y);
% f =@(x,y) -2*sin(x).*cos(y);
% norm_type = 2;
% Nser = [10,50,100,500,1000];
% % Nser = logspace(1,3,5);
% 
% for Nind = 1:length(Nser)
%     N = Nser(Nind);
%     xdom = linspace(0,5,N);
%     ydom = xdom;
%     h = xdom(2)-xdom(1);
% 
%     [Xdom,Ydom] = meshgrid(xdom(2:end-1),ydom(2:end-1));
%     [Xdom_full,Ydom_full] = meshgrid(xdom,ydom);
%     fmat = f(Xdom,Ydom);
%     u_exact = u(Xdom_full,Ydom_full);
% 
%     % Add BCs
%     fmat(1,:) = fmat(1,:) - (1/h^2)*u_exact(1,2:end-1);
%     fmat(end,:) = fmat(end,:) - (1/h^2)*u_exact(end,2:end-1);
%     fmat(:,1) = fmat(:,1) - (1/h^2)*u_exact(2:end-1,1);
%     fmat(:,end) = fmat(:,end) - (1/h^2)*u_exact(2:end-1,end);
% 
%     % Reshape f matrix to vetor
%     fvect = reshape(fmat,(N-2)^2,1);
% 
%     % FDM matrix
%     T = (spdiags(-4*ones(1,N-2),0,N-2,N-2) + spdiags(ones(1,N-3),1,N-2,N-2) + spdiags(ones(1,N-3),-1,N-2,N-2));
%     A1 = kron(speye(N-2),T);
%     A2 = kron(spdiags(ones(1,N-2),1,N-2,N-2),speye(N-2));
%     A3 = kron(spdiags(ones(1,N-2),-1,N-2,N-2),speye(N-2));
%     A = (1/h^2)*(A1 + A2 + A3 );
% 
%     % Solve for unknowns
%     u_appvect = A\fvect;
%     u_approx = reshape(u_appvect,N-2,N-2);
% 
%     % Add BCs
%     u_approx = [u_exact(1,2:end-1); u_approx; u_exact(end,2:end-1)];
%     u_approx = [u_exact(:,1) u_approx u_exact(:,end)];
% 
%     % Calc Error
%     rel_error_2(Nind) = norm(u_exact-u_approx,norm_type)/norm(u_exact,norm_type);
% end
% 
% % Plot error
% figure()
% loglog(Nser,rel_error_2,'.-', 'MarkerSize',10,'Color', '#F67280');
% xlabel('$N$','Interpreter','latex')
% ylabel('Relative Error','Interpreter','latex')
% 
% % Get Slope
% m_2 = polyfit(log(Nser),log(rel_error_2),1)
