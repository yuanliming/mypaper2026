clear; clc;
tic
% parameters
h2_eps = 0.1;
hinf_eps = 0.1;
pole_eps = 0.1;

decay_alpha = 1;

gamma_inf = 3;
tol = 1e-8;

un = 0.1;
up = 1 + un;
lo = 1 - un;

n = 3;
m = 2;
nw = 3;
nz2 = 1;
nzinf = 1;

N = 8;

% nominal system
A0 = [1 4 0;
      1 2 2;
      0 -2 3];

B20 = [1 0;
       1 0;
       0 1];

B1 = eye(3);

C2 = [1 -1 1];
D2 = [1 1];

Cinf = [0 1 -1];
Dinf2 = [1 -1];
Dinf1 = [1 0 0];

% uncertain entries:
% A(1,1), A(2,3), B2(3,2)
A = cell(1,N);
B2 = cell(1,N);

vertex = 0;
for a11_scale = [up lo]
    for a23_scale = [up lo]
        for b32_scale = [up lo]
            vertex = vertex + 1;

            A{vertex} = A0;
            A{vertex}(1,1) = A0(1,1) * a11_scale;
            A{vertex}(2,3) = A0(2,3) * a23_scale;

            B2{vertex} = B20;
            B2{vertex}(3,2) = B20(3,2) * b32_scale;
        end
    end
end

% decision variables
P2 = cell(1,N);
Z2 = cell(1,N);
Pinf = cell(1,N);
P_alpha = cell(1,N);

mat_h2 = cell(1,N);
mat_h2_trace = cell(1,N);
mat_hinf = cell(1,N);
mat_alpha = cell(1,N);
trz = cell(1,N);

for i = 1:N
    P2{i} = sdpvar(n,n,'symmetric');
    Z2{i} = sdpvar(nw,nw,'symmetric');
    Pinf{i} = sdpvar(n,n,'symmetric');
    P_alpha{i} = sdpvar(n,n,'symmetric');
end

% sparse X structure:
% [1 1 0;
%  0 1 0;
%  0 1 1]
X = [sdpvar(1,1) sdpvar(1,1) zeros(1,1);
     zeros(1,1) sdpvar(1,1) zeros(1,1);
     zeros(1,1) sdpvar(1,1) sdpvar(1,1)]; 

% X = [sdpvar(1,1) zeros(1,1) zeros(1,1);
%      zeros(1,1) sdpvar(1,1) zeros(1,1);
%      zeros(1,1) zeros(1,1) sdpvar(1,1)];  % a diagonal X fails to find a feasible solution.


% sparse R structure:
% [1 1 0;
%  0 1 1]
R = [sdpvar(1,1) sdpvar(1,1) zeros(1,1);
     zeros(1,1) sdpvar(1,1) sdpvar(1,1)];

% construct LMIs
for i = 1:N
    AXBR = A{i}*X - B2{i}*R;

    % H2 performance LMI
    Y_h2 = [-X;
             AXBR;
             C2*X - D2*R] ...
             * [h2_eps*eye(n) eye(n) zeros(n,nz2)];

    mat_h2{i} = Y_h2 + Y_h2' ...
              + [zeros(n) P2{i} zeros(n,nz2);
                 P2{i} zeros(n) zeros(n,nz2);
                 zeros(nz2,n) zeros(nz2,n) -eye(nz2)];

    mat_h2_trace{i} = [Z2{i} B1';
                       B1 P2{i}];

    trz{i} = trace(Z2{i});

    % H-infinity constraint with prescribed gamma_inf
    Y_hinf = [-X;
               AXBR;
               Cinf*X - Dinf2*R;
               zeros(nw,n)] ...
               * [hinf_eps*eye(n) eye(n) zeros(n,nzinf) zeros(n,nw)];

    mat_hinf{i} = Y_hinf + Y_hinf' ...
                + [zeros(n) Pinf{i} zeros(n,nzinf) zeros(n,nw);
                   Pinf{i} zeros(n) zeros(n,nzinf) B1;
                   zeros(nzinf,n) zeros(nzinf,n) -gamma_inf*eye(nzinf) Dinf1;
                   zeros(nw,n) B1' Dinf1' -gamma_inf*eye(nw)];

    % pole constraint 1: Re(lambda) < -decay_alpha
    Y_alpha = [-X;
               AXBR;
               zeros(n)] ...
               * [pole_eps*eye(n) eye(n) pole_eps*eye(n)];

    mat_alpha{i} = Y_alpha + Y_alpha' ...
                 + [zeros(n) P_alpha{i} zeros(n);
                    P_alpha{i} zeros(n) P_alpha{i};
                    zeros(n) P_alpha{i} -(1/(2*decay_alpha))*P_alpha{i}];
end

% optimization
gamma_h2 = sdpvar(1);

LMI = [];
for i = 1:N
    LMI = [LMI, ...
           P2{i} >= tol*eye(n), ...
           Z2{i} >= tol*eye(nw), ...
           Pinf{i} >= tol*eye(n), ...
           P_alpha{i} >= tol*eye(n), ...
           mat_h2{i} <= -tol*eye(2*n+nz2), ...
           mat_h2_trace{i} >= tol*eye(n+nw), ...
           mat_hinf{i} <= -tol*eye(2*n+nzinf+nw), ...
           mat_alpha{i} <= -tol*eye(3*n), ...
           trz{i} <= gamma_h2];
end

options = sdpsettings('solver','sdpt3','verbose',1);

sol = optimize(LMI, gamma_h2, options);

if sol.problem ~= 0
    disp(sol.info);
    error('Optimization failed. Try increasing gamma_inf or relaxing pole constraints.');
end

% controller
K_sv = value(R) / value(X)
H2norm_upper = sqrt(value(gamma_h2))
gamma_inf
toc
% verify performance and poles
h2norm_Ksv = zeros(1,N);
hinfnorm_Ksv = zeros(1,N);
cl_poles = zeros(n,N);

for i = 1:N
    Acl = A{i} - B2{i}*K_sv;

    sys_h2 = ss(Acl, B1, C2-D2*K_sv, zeros(nz2,nw));
    sys_hinf = ss(Acl, B1, Cinf-Dinf2*K_sv, Dinf1);

    h2norm_Ksv(i) = norm(sys_h2, 2);
    hinfnorm_Ksv(i) = norm(sys_hinf, inf);
    cl_poles(:,i) = eig(Acl);
end

disp('H2 norms at all vertices:');
disp(h2norm_Ksv);

disp('Maximum H2 norm:');
disp(max(h2norm_Ksv));

disp('H-infinity norms at all vertices:');
disp(hinfnorm_Ksv);

disp('Maximum H-infinity norm:');
disp(max(hinfnorm_Ksv));

disp('Closed-loop poles at all vertices:');
disp(cl_poles);

disp('Maximum real part of closed-loop poles:');
disp(max(real(cl_poles(:))));



cl_poles_h2 = zeros(n,N);

cl_poles_h2(:,1) = [
-0.467803510003123 + 0.593404529050156i;
-0.467803510003123 - 0.593404529050156i;
-5.47983756408749];

cl_poles_h2(:,2) = [
-0.642199266806784 + 0.524428071977687i;
-0.642199266806784 - 0.524428071977687i;
-3.52932960025269];

cl_poles_h2(:,3) = [
-0.436219799791671 + 0.604096391067578i;
-0.436219799791671 - 0.604096391067578i;
-5.54300498451039];

cl_poles_h2(:,4) = [
-0.568818751305473 + 0.564300962588589i;
-0.568818751305473 - 0.564300962588589i;
-3.67609063125531];

cl_poles_h2(:,5) = [
-0.567698116839340 + 0.635080709496452i;
-0.567698116839340 - 0.635080709496452i;
-5.48004835041505];

cl_poles_h2(:,6) = [
-0.741166213865699 + 0.604837839048949i;
-0.741166213865699 - 0.604837839048949i;
-3.53139570613486];

cl_poles_h2(:,7) = [
-0.536138731563082 + 0.639992702553275i;
-0.536138731563082 - 0.639992702553275i;
-5.54316712096757];

cl_poles_h2(:,8) = [
-0.668149931155081 + 0.626758138516824i;
-0.668149931155081 - 0.626758138516824i;
-3.67742827155609];

% plot comparison of poles
figure;
hold on;
box on;
grid on;

% Multi-objective synthesis poles
plot(real(cl_poles(:)), imag(cl_poles(:)), ...
     'b*', 'MarkerSize', 10, 'LineWidth', 1.8);

% H2 synthesis poles
plot(real(cl_poles_h2(:)), imag(cl_poles_h2(:)), ...
     'ro', 'MarkerSize', 10, 'LineWidth', 1.8);

% Pole constraint boundary: Re(lambda) = -decay_alpha
xline(-decay_alpha, 'k--', 'LineWidth', 1.8);

% coordinate labels
xlabel('Real Axis', ...
       'FontSize',14, ...
       'FontName','Times New Roman');

ylabel('Imaginary Axis', ...
       'FontSize',14, ...
       'FontName','Times New Roman');

legend('Multi-objective synthesis', ...
    'Robust $H_2$ synthesis', ...
       'Pole constraint boundary', ...
       'Interpreter','latex', ...
       'Location','best');


% keep equal scaling for pole map
axis equal;

% make sure axes are visible
ax = gca;
ax.FontSize = 14;
ax.FontName = 'Times New Roman';

% optional: enlarge plot range slightly
xlim([-6 1]);
ylim([-2 2]);