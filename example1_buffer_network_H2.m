clear; clc;
tic
%% parameters
h2_eps = 1;
tol = 1e-8;

n = 4;
m = 6;
nw = 4;
nz2 = 10;

%% nominal system
A = [-1 0 0 0;
     0 -2 0 0;
     0 0 -3 0;
     0 0 0 -4];

B2 = [-1 1 0 0 0 0;
       0 -1 -1 1 0 0;
       1 0 1 -1 -1 1;
       0 0 0 0 1 -1];

B1 = eye(4);

C2 = [eye(4);
      zeros(6,4)];

D2 = [zeros(4,6);
      eye(6)];

%% decision variables

P = sdpvar(n,n,'symmetric');
Z = sdpvar(nw,nw,'symmetric');

%% sparse X structure
% 1 0 0 0
% 1 1 0 0
% 1 0 1 0
% 1 1 1 1
X = [sdpvar(1,1) zeros(1,1) zeros(1,1) zeros(1,1);
     sdpvar(1,1) sdpvar(1,1) zeros(1,1) zeros(1,1);
     sdpvar(1,1) zeros(1,1) sdpvar(1,1) zeros(1,1);
     sdpvar(1,1) sdpvar(1,1) sdpvar(1,1) sdpvar(1,1)];

% X = [sdpvar(1,1) 0 0 0;
%      0 sdpvar(1,1) 0 0;
%      0 0 sdpvar(1,1) 0;
%      0 0 0 sdpvar(1,1)];


%% sparse R structure
% [1 0 1 0;
%  1 1 0 0;
%  0 0 0 0;
%  0 0 0 0;
%  0 0 0 0;
%  0 0 0 0]

R = [sdpvar(1,1) zeros(1,1) sdpvar(1,1) zeros(1,1);
     sdpvar(1,1) sdpvar(1,1) zeros(1,1) zeros(1,1);
     zeros(1,4);
     zeros(1,4);
     zeros(1,4);
     zeros(1,4)];

%% H2 synthesis LMI

AXBR = A*X - B2*R;

Y = [-X;
      AXBR;
      C2*X-D2*R] ...
      * [h2_eps*eye(n) eye(n) zeros(n,nz2)];

mat_h2 = Y + Y' ...
       + [zeros(n) P zeros(n,nz2);
          P zeros(n) zeros(n,nz2);
          zeros(nz2,n) zeros(nz2,n) -eye(nz2)];

mat_trace = [Z B1';
             B1 P];

gamma = sdpvar(1);

LMI = [
    P >= tol*eye(n),
    Z >= tol*eye(nw),
    mat_h2 <= -tol*eye(2*n+nz2),
    mat_trace >= tol*eye(n+nw),
    trace(Z) <= gamma
];

%% optimization

options = sdpsettings('solver','sdpt3','verbose',1);

sol = optimize(LMI,gamma,options);

if sol.problem ~= 0
    disp(sol.info);
    error('Optimization failed.');
end

%% controller

K_sv = value(R)/value(X)

H2_upper = sqrt(value(gamma))

toc
%% verification

Acl = A-B2*K_sv;

sys = ss(Acl,B1,C2-D2*K_sv,zeros(nz2,nw));

H2_exact = norm(sys,2)

disp('Closed-loop poles:')
eig(Acl)

disp('Maximum real pole:')
max(real(eig(Acl)))
