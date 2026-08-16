tic
epsilon=0.01;
m=1/400;
d=1/200;
[A1,B21] = MSss(0.3*m, 0.3*d);
[A2,B22] = MSss(0.3*m, -0.3*d);
[A3,B23] = MSss(-0.3*m, 0.3*d);
[A4,B24] = MSss(-0.3*m, -0.3*d);
B1=eye(6);
C=[0 0 0 10^4 0 0;
   0 0 0 0 10^2 0;
   0 0 0 0 0 0;
   0 0 0 0 0 0;];
D=[0;
   0;
   0;
   1];
P1 = sdpvar(6,6,'symmetric');
P2 = sdpvar(6,6,'symmetric');
P3 = sdpvar(6,6,'symmetric');
P4 = sdpvar(6,6,'symmetric');
X = [sdpvar(3,3,'full')   sdpvar(3,3,'full');
    zeros(3,3)     sdpvar(3,3,'full');];

R = [zeros(1,3),sdpvar(1,3,'full')];
gamma = sdpvar(1,1,'full');


LIM1=[-X; A1*X-B21*R; C*X-D*R; zeros(6,6)]*[epsilon*eye(6) eye(6) zeros(6,4) zeros(6,6)] + ([-X; A1*X-B21*R; C*X-D*R; zeros(6,6)]*[epsilon*eye(6) eye(6) zeros(6,4) zeros(6,6)])'...
    +[zeros(6,6) P1 zeros(6,4) zeros(6,6);
       P1 zeros(6,6) zeros(6,4) B1;
       zeros(4,6) zeros(4,6) -gamma*eye(4) zeros(4,6);
       zeros(6,6) B1' zeros(6,4) -gamma*eye(6)];

LIM2=[-X; A2*X-B22*R; C*X-D*R; zeros(6,6)]*[epsilon*eye(6) eye(6) zeros(6,4) zeros(6,6)] + ([-X; A2*X-B22*R; C*X-D*R; zeros(6,6)]*[epsilon*eye(6) eye(6) zeros(6,4) zeros(6,6)])'...
    +[zeros(6,6) P2 zeros(6,4) zeros(6,6);
       P2 zeros(6,6) zeros(6,4) B1;
       zeros(4,6) zeros(4,6) -gamma*eye(4) zeros(4,6);
       zeros(6,6) B1' zeros(6,4) -gamma*eye(6)];

LIM3=[-X; A3*X-B23*R; C*X-D*R; zeros(6,6)]*[epsilon*eye(6) eye(6) zeros(6,4) zeros(6,6)] + ([-X; A3*X-B23*R; C*X-D*R; zeros(6,6)]*[epsilon*eye(6) eye(6) zeros(6,4) zeros(6,6)])'...
    +[zeros(6,6) P3 zeros(6,4) zeros(6,6);
       P3 zeros(6,6) zeros(6,4) B1;
       zeros(4,6) zeros(4,6) -gamma*eye(4) zeros(4,6);
       zeros(6,6) B1' zeros(6,4) -gamma*eye(6)];

LIM4=[-X; A4*X-B24*R; C*X-D*R; zeros(6,6)]*[epsilon*eye(6) eye(6) zeros(6,4) zeros(6,6)] + ([-X; A4*X-B24*R; C*X-D*R; zeros(6,6)]*[epsilon*eye(6) eye(6) zeros(6,4) zeros(6,6)])'...
    +[zeros(6,6) P4 zeros(6,4) zeros(6,6);
       P4 zeros(6,6) zeros(6,4) B1;
       zeros(4,6) zeros(4,6) -gamma*eye(4) zeros(4,6);
       zeros(6,6) B1' zeros(6,4) -gamma*eye(6)];

problem=[LIM1<=0,LIM2<=0,LIM3<=0,LIM4<=0,P1>=0,P2>=0,P3>=0,P4>=0,gamma>=0];

options=sdpsettings('solver','sdpt3','verbose',0);
optimize(problem, gamma,options); 

K_svinf = value(R)*inv(value(X)) 
gamma=(value(gamma))
toc
vertex1=norm(ss(A1-B21*K_svinf,B1,C-D*K_svinf,zeros(4,6)),inf)
vertex2=norm(ss(A2-B22*K_svinf,B1,C-D*K_svinf,zeros(4,6)),inf)
vertex3=norm(ss(A3-B23*K_svinf,B1,C-D*K_svinf,zeros(4,6)),inf)
vertex4=norm(ss(A4-B24*K_svinf,B1,C-D*K_svinf,zeros(4,6)),inf)

function [A,B2] = MSss(dm, dd)
m=1/400;
d=1/200;
a1=-125;a2=-75;a3=-15;

A=[0 1 0 0 0 0;
   0 0 1 0 0 0;
   a1 a2 a3 0 0 0;
   0 0 0 0 1 0;
   0 0 0 0 0 1;
   a1*dm/(m+dm) a2*dm/(m+dm) (a3*dm+dd)/(m+dm) 0 0 -(d+dd)/(m+dm)];
B2=[0;0;0;0;0;-1/(m+dm)];
end
