clear all;
clc;
N=20;
K=1*[-0.56   0       0.56    0;
          0  -0.56    0       0.56;
      0.218   0.218   0.218   0.218];
a=0.3490000;

V=[1;1;1];

e = V;
e=e/norm(e);
cf = @(x)confun(x,e,K);
of = @(x)objfun(x,K);
% A=[];b=[];Aeq=[];beq=[];
A=-e'*K;b=0;Aeq=[0 -e(3) e(2);e(3) 0 -e(1);-e(2) e(1) 0]*K;beq=[0;0;0];
lb=[-a,-a,-a,-a];ub=[a,a,a,a];
options = optimoptions(@ fmincon,'Algorithm','interior-point','Display','iter');
% x0=[0;0;0;0]; %对解决方案进行初步猜测
c1=-e(1)+0.5*e(3);c2=-e(2)+0.5*e(3);c3= e(1)+0.5*e(3);c4= e(2)+0.5*e(3);
x0=0.1*[c1;c2;c3;c4];
[x,fval,exitflag,output,lambda,grad,hessian] = fmincon(of,x0,A,b,Aeq,beq,lb,ub,[],options);
U=K*x;

% load seamount
figure,
plot3(U(1,:),U(2,:),U(3,:),'r.');hold on;

figure,
plot3(U(1,:),U(2,:),U(3,:),'r.');hold on;
plot3(V(1),V(2),V(3),'*');
% surf(X,Y,Z);