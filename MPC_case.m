%% 加载 optim package,若使用matlab，则注释掉此行

pkg load optim;

%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 第一步，定义状态空间矩阵

%% 定义状态矩阵 A, n x n 矩阵

A = [1 0.1; -1 2];

n= size (A,1);

%% 定义输入矩阵 B, n x p 矩阵

B = [ 0.2 1; 0.5 2];

p = size(B,2);

%% 定义Q矩阵，n x n 矩阵

Q=[100 0;0 1];

%% 定义F矩阵，n x n 矩阵

F=[100 0;0 1];

%% 定义R矩阵，p x p 矩阵

R=[1 0 ;0 .1];

%% 定义step数量k

k_steps=100; 

%% 定义矩阵 X_K， n x k 矩 阵

X_K = zeros(n,k_steps);

%% 初始状态变量值， n x 1 向量

X_K(:,1) =[20;-20];

%% 定义输入矩阵 U_K， p x k 矩阵

U_K=zeros(p,k_steps);


%% 定义预测区间K

N=5;

%% Call MPC_Matrices 函数 求得 E,H矩阵 

[E,H]=MPC_Matrices(A,B,Q,R,F,N);


%% 计算每一步的状态变量的值

for k = 1 : k_steps 

%% 求得U_K(:,k)

U_K(:,k) = Prediction(X_K(:,k),E,H,N,p);

%% 计算第k+1步时状态变量的值

X_K(:,k+1)=(A*X_K(:,k)+B*U_K(:,k));

end


%% 绘制状态变量和输入的变化

subplot  (2, 1, 1);

hold;

for i =1 :size (X_K,1)

plot (X_K(i,:));

end

legend("x1","x2")

hold off;

subplot (2, 1, 2);

hold;

for i =1 : size (U_K,1)

plot (U_K(i,:));

end

legend("u1","u2")

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

~~~~~~~~~~~~~~~~~~~~~~~MPC_Matrices.m~~~~~~~~~~~~~~~~~~~~

function  [E , H]=MPC_Matrices(A,B,Q,R,F,N)


n=size(A,1);   % A 是 n x n 矩阵, 得到 n

p=size(B,2);   % B 是 n x p 矩阵, 得到 p

%%%%%%%%%%%%

M=[eye(n);zeros(N*n,n)]; % 初始化 M 矩阵. M 矩阵是 (N+1)n x n的， 

                         % 它上面是 n x n 个 "I", 这一步先把下半部

                         % 分写成 0 

C=zeros((N+1)*n,N*p); % 初始化 C 矩阵, 这一步令它有 (N+1)n x NP 个 0

% 定义M 和 C 

tmp=eye(n);  %定义一个n x n 的 I 矩阵

%　更新Ｍ和C

for i=1:N % 循环，i 从 1到 N

    rows =i*n+(1:n); %定义当前行数，从i x n开始，共n行 

    C(rows,:)=[tmp*B,C(rows-n, 1:end-p)]; %将c矩阵填满

    tmp= A*tmp; %每一次将tmp左乘一次A

    M(rows,:)=tmp; %将M矩阵写满

end 


% 定义Q_bar和R_bar

Q_bar = kron(eye(N),Q);

Q_bar = blkdiag(Q_bar,F);

R_bar = kron(eye(N),R); 


% 计算G, E, H

G=M'*Q_bar*M; % G: n x n

E=C'*Q_bar*M; % E: NP x n

H=C'*Q_bar*C+R_bar; % NP x NP 


end

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

~~~~~~~~~~~~~~~~~~~~~~~Prediction.m~~~~~~~~~~~~~~~~~~~~~

function u_k= Prediction(x_k,E,H,N,p)

U_k = zeros(N*p,1); % NP x 1

U_k = quadprog(H,E*x_k);

u_k = U_k(1:p,1); % 取第一个结果

end

