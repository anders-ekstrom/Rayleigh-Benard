%% 
%       2D incompressible Rayleigh-Benard Convection, Boussinesq Approx.
%       second order space and time (second order central spatial, RK2
%       temporal, pressure projection)
%       periodic in x, no slip at y=(0,H), theta dirichlet BC at y=(0,H)
%       u is nx x ny, stored at x faces and y centers
%       v is nx x (ny+1), stored at y faces
%       theta is nx x (ny+2), stored at cell centers with ghost points
%%
clc
clear
clf
close all

%domain size
L = 2;
H = 0.5;

%cells in x and y
nx = 400;
ny = 200;

%grid
dx = L/nx;
dy = H/ny;

%for plotting
x_center = ((1:nx)-1/2).*dx;
y_center = ((1:ny)-1/2).*dy;

%flow paramaters
Pr = 10;
Ra = 1e7;

%simulation parameters
t_end = 20; 
cfl = 0.4;
plot_dt = 0.01;
next_frame = 0;

%boundary condtions
theta_bottom = 10;
theta_top = -10;

%initial condtions
u = zeros(nx,ny);
v = zeros(nx,ny+1);
theta = zeros(nx,ny+2);
theta(:,2:ny+1) = ones(nx,ny) .* 0.5 * (theta_bottom + theta_top);
theta = theta_bc(theta,theta_bottom,theta_top);

%see functions below, creating x and y operator matrices
[ddx_forward,ddx_backward,d2dx2] = x_operator(nx,dx);
[ddy_ctf,ddy_ftc,d2dy2_theta,d2dy2_v,d2dy2_u,d2dy2_poisson] = y_operator(ny,dy);

% make poisson operator for projection
Ix = speye(nx);
Iy = speye(ny);
A = kron(Iy,d2dx2) + kron(d2dy2_poisson,Ix);
A(1,:) = 0;
A(1,1) = 1;
A = decomposition(A,'lu');

%main time loop 
t = 0;
iter = 0;

while t < t_end
    iter = iter + 1;

    %set stable dt (advection, diffusion, and simulation constraints)
    dt = min([cfl*dx / (max(abs(u(:)))+eps), cfl*dy / (max(abs(v(:)))+eps),...
                    1/4*min(dx^2,dy^2)/(max(sqrt(Pr/Ra),sqrt(1/(Ra*Pr)))), t_end - t]);

    %predictor step (euler)
    [fu,fv,ftheta] = RHS(u,v,theta,Pr,Ra,nx,ny,dx,dy,theta_bottom,theta_top,d2dx2,d2dy2_u,d2dy2_v,d2dy2_theta);
    u_euler = u + dt*fu;
    v_euler = v + dt*fv;
    theta_euler = theta;
    theta_euler(:,2:ny+1) = theta(:,2:ny+1) + dt*ftheta;
    %apply bc
    theta_euler = theta_bc(theta_euler,theta_bottom,theta_top);
    v_euler(:,1) = 0;
    v_euler(:,end) = 0;

    %projection
    [u_proj_euler, v_proj_euler] = velocity_project(u_euler,v_euler,nx,ny,dt,A,ddx_backward,ddx_forward,ddy_ftc,ddy_ctf);

    %boundary condtions again
    v_proj_euler(:,1) = 0;
    v_proj_euler(:,end) = 0;

    %corrector step
    [cu, cv, ctheta] = RHS(u_proj_euler,v_proj_euler,theta_euler,Pr,Ra,nx,ny,dx,dy,theta_bottom,theta_top,d2dx2,d2dy2_u,d2dy2_v,d2dy2_theta);
    u_np1 = u + 0.5*dt*(fu + cu);
    v_np1 = v + 0.5*dt*(fv + cv);
    theta_np1 = theta;
    theta_np1(:,2:ny+1) = theta(:,2:ny+1) + 0.5*dt*(ftheta + ctheta);
    
    %projection again
    [u_np1, v_np1] = velocity_project(u_np1, v_np1, nx, ny, dt, A, ddx_backward, ddx_forward, ddy_ftc, ddy_ctf);

    %boundary conditions yet again
    theta_np1 = theta_bc(theta_np1,theta_bottom,theta_top);
    v_np1(:,1) = 0;
    v_np1(:,end) = 0;

    u = u_np1;
    v = v_np1;
    theta = theta_np1;
    t = t + dt;

    if t >= next_frame
        next_frame = next_frame + plot_dt;
        theta_plot = theta(:,2:ny+1);
        pcolor(x_center,y_center,theta_plot')
        shading interp
        axis xy
        axis equal
        xlim([0 L])
        ylim([0 H])
        colorbar
        clim([theta_top theta_bottom])
        colormap(cool_hot())
        title(sprintf('\\theta at t = %.2f ', t))
        drawnow
    end
end


%% Functions
function [ddx_forward,ddx_backward,d2dx2] = x_operator(nx,dx)
    % forward, backward, central, and laplace matrix operators for periodic x
    % Spdiags documentation
    % S = spdiags(Bin,d,m,n) creates an m-by-n sparse matrix S by taking the columns of Bin and placing them along the diagonals specified by d.
    
    shift_forward = spdiags(ones(nx,1),1,nx,nx);
    shift_forward(nx,1) = 1;
    
    shift_backward = spdiags(ones(nx,1),-1,nx,nx);
    shift_backward(1,nx) = 1;
    
    I = speye(nx,nx);
    
    ddx_forward = (shift_forward - I)./dx;
    ddx_backward = (I - shift_backward)./dx;
    d2dx2 = (shift_forward - 2.*I + shift_backward)./(dx^2);
end

function [ddy_ctf,ddy_ftc,d2dy2_theta,d2dy2_v,d2dy2_u,d2dy2_poisson] = y_operator(ny,dy)
    % center to face
    ddy_ctf = spalloc(ny+1,ny,2*ny-2);
    for i = 2:ny
        ddy_ctf(i,i) = 1/dy;
        ddy_ctf(i,i-1) = -1/dy;
    end

    % face to center
    ddy_ftc = spalloc(ny,ny+1,2*ny);
    for i= 1:ny
        ddy_ftc(i,i) = -1/dy;
        ddy_ftc(i,i+1) = 1/dy;
    end
    
    % second derivitve interior theta. shifted cause ghost cells
    d2dy2_theta = spalloc(ny,ny+2,3*ny);
    for i = 1:ny
        d2dy2_theta(i,i) = 1/dy^2;
        d2dy2_theta(i,i+1) = -2/dy^2;
        d2dy2_theta(i,i+2) = 1/dy^2;
    end

    % standard operator, only set interior as bc at wall is v=0
    d2dy2_v = spalloc(ny+1,ny+1,ny*3+1);
    for i = 2:ny
        d2dy2_v(i,i-1) = 1/dy^2;
        d2dy2_v(i,i) = -2/dy^2;
        d2dy2_v(i,i+1) = 1/dy^2;
    end

    % in u momentum eqn. modified stencil to enforce u=0 at wall
    d2dy2_u = spalloc(ny,ny,ny*3);
    for i = 2:ny-1
        d2dy2_u(i,i-1) = 1/dy^2;
        d2dy2_u(i,i) = -2/dy^2;
        d2dy2_u(i,i+1) = 1/dy^2;
    end
    d2dy2_u(1,1) = -3/dy^2;
    d2dy2_u(1,2) = 1/dy^2;
    d2dy2_u(ny,ny) = -3/dy^2;
    d2dy2_u(ny,ny-1) = 1/dy^2;

    % for poisson equation. modeified stencil for nueman bc.
    d2dy2_poisson = spalloc(ny,ny,ny*3);
    for i = 2:ny-1
        d2dy2_poisson(i,i-1) = 1/dy^2;
        d2dy2_poisson(i,i) = -2/dy^2;
        d2dy2_poisson(i,i+1) = 1/dy^2;
    end
    d2dy2_poisson(1,1) = -1/dy^2;
    d2dy2_poisson(1,2) = 1/dy^2;
    d2dy2_poisson(ny,ny) = -1/dy^2;
    d2dy2_poisson(ny,ny-1) = 1/dy^2;
end

function theta = theta_bc(theta,theta_bottom,theta_top)
    % theta boundary condtions using ghost cells
    theta(:,1) = 2*theta_bottom - theta(:,2);
    theta(:,end) = 2*theta_top - theta(:,end-1);
end



function [RHS_u, RHS_v, RHS_theta] = RHS(u,v,theta,Pr,Ra,nx,ny,dx,dy,theta_bottom,theta_top,d2dx2,d2dy2_u,d2dy2_v,d2dy2_theta)
    theta_interior = theta(:,2:ny+1);

    % more operators
    x_forward = @(Q) circshift(Q,[-1,0]);
    x_backward = @(Q) circshift(Q,[1,0]);
    x_avg = @(Q) 1/2 * (Q+x_forward(Q));
    ddx_flux = @(Q) (Q-x_backward(Q))/dx;
    ddy_flux = @(Q) (Q(:,2:end)-Q(:,1:end-1))/dy;


    %since staggared grid need to evaluate u and v at same places for some
    %terms
    u_yface = zeros(nx,ny+1);
    u_yface(:,2:ny) = 1/2 * (u(:,1:ny-1) + u(:,2:ny));
    %enforce bc
    u_yface(:,1) = 0;
    u_yface(:,end) = 0;

    v_xmid = x_avg(v);

    uv_corner = u_yface.* v_xmid;

    % laplacians / diffusion
    laplacian_u = d2dx2*u + u*d2dy2_u';
    laplacian_v = d2dx2*v + v*d2dy2_v';
    laplacian_theta = d2dx2*theta_interior + theta*d2dy2_theta';

    % u momentum equation
    % du/dt = -d(uu)/dx - d(uv)/dy) + (pr/ra)^0.5 * (d2(u)/dx2 + d2(u)/dy2)

    u_xmid = x_avg(u);
    duudx = ddx_flux(u_xmid.^2);

    %make sure evaluted on x face (where u is )
    duvdy = ddy_flux(uv_corner(:,1:ny+1));

    RHS_u = - duudx - duvdy + (Pr/Ra)^0.5 * laplacian_u;

    % v momentum equation
    % dv/dt = -d(uv)/dx - d(vv)/dy + (pr/ra)^0.5 * (d2(v)/dx2 + d2(v)/dy2)
    % + theta

    %make sure is evaluated on y face (where v is)
    duvdx = ddx_flux(uv_corner);
    
    v_ymid = 1/2 * (v(:,1:ny) + v(:,2:ny+1));
    vv_ymid = v_ymid.^2;
    % maintains no slip bc
    dvvdy = zeros(nx, ny+1);
    dvvdy(:,2:ny) = (vv_ymid(:,2:ny) - vv_ymid(:,1:ny-1))/dy;

    theta_yface = zeros(nx,ny+1);
    theta_yface(:,2:ny) = 1/2 * (theta_interior(:,1:ny-1) + theta_interior(:,2:ny));

    RHS_v = -duvdx - dvvdy + (Pr/Ra)^0.5 * laplacian_v + theta_yface;
    %enforce bc
    RHS_v(:,1) = 0;
    RHS_v(:,end) = 0;


    % theta equation: dtheta/dt = -d(u theta)/dx - d(v theta)/dy +
    % (1/(Pr*Ra)^0.5 * (d2(theta)/dx2 + d2(theta)/dy2)

    theta_xmid = x_avg(theta_interior);
    utheta = u.*theta_xmid;
    duthetadx = ddx_flux(utheta);

    theta_yface = zeros(nx,ny+1);
    theta_yface(:,1) = theta_bottom;
    theta_yface(:,end) = theta_top;
    theta_yface(:,2:ny) = 1/2 * (theta_interior(:,1:ny-1) + theta_interior(:,2:ny));

    vtheta = v.* theta_yface;
    dvthetady = ddy_flux(vtheta);

    RHS_theta = -duthetadx - dvthetady + (1/(Pr*Ra))^0.5 * laplacian_theta;
end
    
    
function [u_proj, v_proj] = velocity_project(u_star,v_star,nx,ny,dt,A,ddx_backward,ddx_forward,ddy_ftc,ddy_ctf)
    divergence = ddx_backward*u_star + v_star*ddy_ftc';
    RHS = 1/dt*divergence;
    RHS_vector = RHS(:);
    RHS_vector(1) = 0;

    phi_vector = A \ RHS_vector;
    phi = reshape(phi_vector,nx,ny);

    dphidx = ddx_forward *phi;
    dphidy = phi * ddy_ctf';

    u_proj = u_star - dt*dphidx;
    v_proj = v_star - dt*dphidy;
end

function cmap = cool_hot()
    num = 512;
    %blue
    cold = [0.10 0.30 0.95];
    %white
    mid  = [1.00 1.00 1.00];
    %red
    hot  = [0.95 0.25 0.05];

    c1 = [linspace(cold(1),mid(1),num/2)', linspace(cold(2),mid(2),num/2)', linspace(cold(3),mid(3),num/2)'];
    c2 = [linspace(mid(1),hot(1),num/2)',  linspace(mid(2),hot(2),num/2)',  linspace(mid(3),hot(3),num/2)'];

    cmap = [c1; c2];
end

