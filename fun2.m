function flag=fun2(varargin)
%% ��������뺯�����ж�ֱ���Ƿ���ĳ��ƽ����
% ����Ϊ���㵽�߽������������Ϊ����ֵ

% �������ε���������ֱ�ΪABC����Ҫ�жϵĵ�ΪP�������м��ַ����������жϣ�
% 
% ����1����������жϡ�
% ��P�����������ߵ���ɸ������������Σ�PAB��PAC��PBC��
% �ж����ǵ����֮���Ƿ���ABC�������ȡ�����������Ҫ�жϵĵ��������������ڲ���
% ������������������ⲿ��
% 
% ����2���ýǶ����жϡ�
% ��APB��BPC��CPA֮���Ƿ����360�ȣ������������Ҫ�жϵĵ��������������ڲ���
% ������������������ⲿ��
% 
% ����3�������������жϡ�
% �� ��ABC ������߰�һ�������ߣ�˳ʱ�����ʱ�룩��
% �жϵ� P �Ƿ��ڸñߵ�ĳ�ࣨ�Ҳ����ࣩ��
% ���� P �������ߵ�ͬ�࣬��� P �� ��ABC �ڡ�
% �� ��ABC ������߰�һ�������ߣ�˳ʱ�����ʱ�룩��
% �жϵ� P �Ƿ��ڸñߵ�ĳ�ࣨ�Ҳ����ࣩ������ P �������ߵ�ͬ�࣬��� P �� ��ABC �ڡ� 
%%
if nargin==3
    ang=fun4(varargin{1},varargin{2})+fun4(varargin{2},varargin{3})+fun4(varargin{3},varargin{1})
    if ang>=6.2
        flag=true;
    else
        flag=false;
    end
elseif nargin==4
    ang=fun4(varargin{1},varargin{2})+fun4(varargin{2},varargin{3})+fun4(varargin{3},varargin{4})+fun4(varargin{4},varargin{1})
    if ang>=6.2
        flag=true;
    else
        flag=false;
    end
end
% nargin
% varargin{:}
% varargout=varargin;
% nargout
% varargout{:}
end