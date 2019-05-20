function flag=fun2(varargin)
%% 变参数输入函数，判断直线是否在某个平面内
% 输入为交点到边界点的向量，输出为布尔值

% 设三角形的三个顶点分别为ABC，需要判断的点为P。这里有几种方法来进行判断：
% 
% 方法1：用面积来判断。
% 点P与三角形三边的组成个了三个三角形：PAB，PAC，PBC。
% 判断它们的面积之和是否与ABC的面积相等。如果相等则需要判断的点在三角形区域内部，
% 否则就在三角形区域外部。
% 
% 方法2：用角度来判断。
% 角APB，BPC，CPA之和是否等于360度，如果等于则需要判断的点在三角形区域内部，
% 否则就在三角形区域外部。
% 
% 方法3：用向量积来判断。
% 沿 △ABC 各有向边按一定方向走（顺时针或逆时针），
% 判断点 P 是否在该边的某侧（右侧或左侧），
% 若点 P 在三条边的同侧，则点 P 在 △ABC 内。
% 沿 △ABC 各有向边按一定方向走（顺时针或逆时针），
% 判断点 P 是否在该边的某侧（右侧或左侧），若点 P 在三条边的同侧，则点 P 在 △ABC 内。 
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