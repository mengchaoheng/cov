function M = fun1(A,B,C,D,u)
%% ��ֱ����ƽ��Ľ�������
% ��֪�ռ�ֱ��L��(x-a)/m=(x-b)/n=(z-c)/p�Ϳռ�ƽ��У�Ax+By+Cz+D=0;
% ��ֱ��L��ƽ��еĽ�������ꡣ
% ��ֱ�߷��̸�д�ɲ�����ʽ����(x-a)/m=(x-b)/n=(z-c)/p=t��
% ��x=mt+a��y=nt+b��z=pt+c������ƽ��еķ��̵ã�
% A(mt+a)+B(nt+b)+C(pt+c)+D=0
% �ɴ˽��t=-(Aa+Bb+Cc+D)/(Am+Bn+Cp)
% �ٴ���������̼��ý��������(x��y��z).
%%
m=u(1);
n=u(2);
p=u(3);
a=0;
b=0;
c=0;
t=-(A*a+B*b+C*c+D)/(A*m+B*n+C*p);
M=[t*m;t*n;t*p];
end