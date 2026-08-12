clear; clc;

% parameters
h2_eps = 0.1;
tol = 1e-8;

un = 0.1;
up = 1 + un;
lo = 1 - un;

n = 3;
m = 2;
nw = 3;
nz2 = 1;

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

gamma_inf = 3;
decay_alpha = 0.5;
disk_radius = 5;

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
mat_h2 = cell(1,N);
mat_h2_trace = cell(1,N);
trz = cell(1,N);

for i = 1:N
    P2{i} = sdpvar(n,n,'symmetric');
    Z2{i} = sdpvar(nw,nw,'symmetric');
end

% sparse X structure:
% [1 1 0;
%  0 1 0;
%  0 1 1]
X = [sdpvar(1,1) sdpvar(1,1) zeros(1,1);
     zeros(1,1) sdpvar(1,1) zeros(1,1);
     zeros(1,1) sdpvar(1,1) sdpvar(1,1)];

% sparse R structure:
% [1 1 0;
%  0 1 1]
R = [sdpvar(1,1) sdpvar(1,1) zeros(1,1);
     zeros(1,1) sdpvar(1,1) sdpvar(1,1)];

% construct H2 LMIs
for i = 1:N
    AXBR = A{i}*X - B2{i}*R;

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
end

% optimization
gamma_h2 = sdpvar(1);

LMI = [];
for i = 1:N
    LMI = [LMI, ...
           P2{i} >= tol*eye(n), ...
           Z2{i} >= tol*eye(nw), ...
           mat_h2{i} <= -tol*eye(2*n+nz2), ...
           mat_h2_trace{i} >= tol*eye(n+nw), ...
           trz{i} <= gamma_h2];
end

options = sdpsettings('solver','sdpt3','verbose',1);

sol = optimize(LMI, gamma_h2, options);

if sol.problem ~= 0
    disp(sol.info);
    error('Optimization failed.');
end

% controller
K_sv = value(R) / value(X)
H2norm_upper = sqrt(value(gamma_h2))

K_sv = value(R) / value(X)
H2norm_upper = sqrt(value(gamma_h2))
gamma_inf

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


