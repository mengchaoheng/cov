clear all;
clc;
% close all;
N=100;
K=1*[-0.56   0       0.56    0;
          0  -0.56    0       0.56;
      0.218   0.218   0.218   0.218];
a=0.3490000;
% 平面方程
A1=0.3893;B1=0.3893;C1=1;D1=-0.3043;
A2=1;B2=0;C2=0;D2=-0.39088;
A3=0;B3=1;C3=0;D3=-0.39088;
% 平面的点
a1=[0;0;0.3043];b1=[0.39088;0;0.1521];c1=[0.39088;0.3088;0];d1=[0;0.39088;0.1521];
a2=[0.39088;0;0];b2=[0.39088;0;0.1521];c2=[0.39088;0.39088;0];%
a3=[0;0.39088;0];b3=[0;0.39088;0.1521];c3=[0.39088;0.39088;0];

[X,Y,Z] = sphere(N);
U=zeros((N+1)^2,3);
fval=zeros((N+1)^2,1);

for i=1:(N+1)^2
e = [X(i);Y(i);Z(i)];
e1=fun3(e);
flag1=0;
flag2=0;
flag3=0;
M1 = fun1(A1,B1,C1,D1,e1);
if (M1(1)<=0.39088 && M1(1)>=0 && M1(2)<=0.39088 && M1(2)>=0 && M1(3)<=0.3043 && M1(3)>=0)
    flag1=1;
end

if(dot([A2;B2;C2],e1)>1e-5)
    M2 = fun1(A2,B2,C2,D2,e1);
    if (0.39088*M2(3)+0.1521*M2(2)-0.39088*0.1521<0 && M2(2)<0.39088 && M2(2)>=0 && M2(3)<0.1521 && M2(3)>=0)
        flag2=1;
    end
end

if(dot([A3;B3;C3],e1)>1e-5)
    M3 = fun1(A3,B3,C3,D3,e1);
    if (M3(1)<0.39088 && M3(1)>=0 && 0.39088*M3(3)+0.1521*M3(1)-0.39088*0.1521<0 && M3(3)<0.1521 && M3(3)>=0)
        flag3=1;
    end
end

if flag1
    [~,~,r_m] = cart2sph(M1(1),M1(2),M1(3));
    [azimuth,elevation,r_e] = cart2sph(e(1),e(2),e(3));
    if r_e>r_m
        r_e=r_m;
    end
    [U(i,1),U(i,2),U(i,3)]=sph2cart(azimuth,elevation,r_e);
end
if flag2
    [~,~,r_m] = cart2sph(M2(1),M2(2),M2(3));
    [azimuth,elevation,r_e] = cart2sph(e(1),e(2),e(3));
    if r_e>r_m
        r_e=r_m;
    end
    [U(i,1),U(i,2),U(i,3)]=sph2cart(azimuth,elevation,r_e);
end
if flag3
    [~,~,r_m] = cart2sph(M3(1),M3(2),M3(3));
    [azimuth,elevation,r_e] = cart2sph(e(1),e(2),e(3));
    if r_e>r_m
        r_e=r_m;
    end
    [U(i,1),U(i,2),U(i,3)]=sph2cart(azimuth,elevation,r_e);
    fval(i)=r_e;
end
end
U1=zeros(N+1,N+1);
U2=zeros(N+1,N+1);
U3=zeros(N+1,N+1);
for i=1:(N+1)^2
    U1(i)=U(i,1);
    U2(i)=U(i,2);
    U3(i)=U(i,3);
end
figure,
plot3(U(:,1),U(:,2),U(:,3),'r.');hold on;
C=sqrt(U1.*U1+U2.*U2+U3.*U3);
figure,
surf(U1,U2,U3,C);% 'FaceColor','b','FaceAlpha',0.5,,'EdgeColor','none'
% hold on;
% plot3(0,0,0.3043,'r.');hold on;
% plot3(0.3088,0,0.1521,'r.');hold on;
% plot3(0.3088,0.3088,0,'r.');hold on;
% plot3(0,0.3088,0.1521,'r.');hold on;
% plot3(0.3088,0,0,'r.');hold on;
% plot3(0,0.3088,0,'r.');hold on;
% plot3(0,0.3088,0.1521,'r.');
% x1=linspace(0,1,10);
% y1=linspace(0,1,10);
% [X1,Y1]=meshgrid(x1,y1);
% Z1=0.3043 - X1.*0.3893 - Y1.*0.3893;
% mesh(X1,Y1,Z1)
% y1=linspace(0,1,10);
% z1=linspace(0,1,10);
% [Y1,Z1]=meshgrid(y1,z1);
% X1=ones(10).*0.39088;
% mesh(X1,Y1,Z1)
% x1=linspace(0,1,10);
% z1=linspace(0,1,10);
% [X1,Z1]=meshgrid(x1,z1);
% Y1=ones(10).*0.39088;
% mesh(X1,Y1,Z1)

colorbar
C1=sqrt(X.*X+Y.*Y+Z.*Z);
figure,
surf(X,Y,Z,C1);
colorbar
