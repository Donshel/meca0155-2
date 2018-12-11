%% Beam

beam = struct;

% Properties
beam.dim = [1.25; 0.02; 0.08]; % m
beam.vol = prod(beam.dim); % m^3
beam.density = 2690; % kg / m^3
beam.mass = beam.density * beam.vol; % kg

beam.E = 6.45 * 10^10; % N / m^
beam.v = 0.39;

% Position
beam.pos = [0; 0; 0];

% Inertia moment
beam.J = zeros(3, 1);
for i = 1:3
    temp = (1:3)'; temp(i) = [];
    beam.J(i) = beam.mass * sum(beam.dim(temp).^2) / 12;
end

%% Cylinder

cyl = struct;

% Properties
cyl.dim = [0.051; 0.061; 0.051]; % m
cyl.vol = prod(cyl.dim) * pi / 4; % m^3
cyl.mass = 1; % kg
cyl.density = cyl.mass / cyl.vol; % kg / m^3

% Position
cyl.pos = [0.3; (cyl.dim(2) + beam.dim(2)) / 2; 0]; % m

% Inertia moment
cyl.J = zeros(3, 1);
for i = [1 3]
    cyl.J(i) = cyl.mass * (cyl.dim(2)^2 / 12 + sum((cyl.dim([1;3]) / 2).^2) / 8 );
end
temp = [1; 3]; % i = 2
cyl.J(2) = cyl.mass * sum((cyl.dim(temp) / 2).^2) / 4;

%% CM

cm = (beam.pos * beam.mass + cyl.pos * cyl.mass) / (beam.mass + cyl.mass);

% Transports
d = beam.pos - cm;
beam.trsp = zeros(3, 1);
for i = 1:3
    temp = (1:3)';
    temp(i) = [];
    beam.trsp(i) = beam.mass * sum(d(temp).^2);
end

d = cyl.pos - cm;
cyl.trsp = zeros(3, 1);
for i = 1:3
    temp = (1:3)';
    temp(i) = [];
    cyl.trsp(i) = cyl.mass * sum(d(temp).^2);
end

%% Supports

string = struct;

% Properties
string(1).stiff = 10500; % N / m
string(2).stiff = string(1).stiff; % N / m

% Position
string(1).pos = [-0.5; - beam.dim(2) / 2; 0]; % m
string(2).pos = [0.4; - beam.dim(2) / 2; 0]; % m

%% Symbolic computation

% Parameters
syms m_eq J_eq x1 x2 k1 k2;

% Generelized coordinates
syms t y(t) theta(t);

f = sym('f', [2 3]);
f(1, 1) = y;
f(2, 1) = theta;
for i = 1:size(f, 1)
    for j = 2:size(f, 2)
        f(i, j) = diff(f(i, j - 1), t);
    end
end

q = sym('q', size(f));
h1 = q(1, 1) + x1 * q(2, 1); % tan(theta) ~ theta
h2 = q(1, 1) + x2 * q(2, 1);

% Lagrangian
V = k1 * h1^2 / 2 + k2 * h2^2 / 2;
T = m_eq * q(1, 2)^2 / 2 + J_eq * q(2, 2)^2 / 2;
L = T - V;

% Matrices
mat = cell(3, 1);
for i = 1:3
    mat{i} = sym('m', [size(f, 1) size(f, 1)]);
end

E = sym([size(f, 1) 1]);
for i = 1:size(f, 1)
    E(i) = diff(subs(diff(L, q(i, 2)), q, f), t) - subs(diff(L, q(i, 1)), q, f);
    for j = 1:size(f, 1)
        for k = 1:3
            mat{k}(i, j) = 0;
            [c, t] = coeffs(E(i), f(j, k));
            for l = 1:length(t)
                if isequaln(t(l), f(j,k))
                    mat{k}(i, j) = c(l);
                    break;
                end
            end
        end
    end
end

K = mat{1};
M = mat{3};

%% Substitution

m_eq = beam.mass + cyl.mass;
J_eq = beam.J(3) + beam.trsp(3) + cyl.J(3) + cyl.trsp(3);
x1 = string(1).pos(1) - cm(1);
x2 = string(2).pos(1) - cm(1);
k1 = string(1).stiff;
k2 = string(2).stiff;

subst = struct;

subst.K = double(subs(K));
subst.M = double(subs(M));

%% Eigenvalues and eigenvectors

[subst.V, subst.D2] = eig(subst.K, subst.M);
subst.D = (subst.D2).^(1/2);

for i = 1:size(subst.V, 1)
    subst.D(i, i) = subst.D(i, i) / (2 * pi);
    subst.V(:, i) = subst.V(:, i) / subst.V(i, i);
end

subst.f = diag(subst.D);

%% Node

subst.x = -subst.V(1, 2) / subst.V(2, 2) + cm(1);

%% Compare

freq = load('resources/txt/Project_2018_freq.txt');
err = subst.f(1:2) - freq(1:2, 1);
rel_err = err ./ freq(1:2, 1);

%% Display

disp('K =');
disp(subst.K);
disp('M =');
disp(subst.M);
disp('freq =');
disp(subst.f);
disp('modes =');
disp(subst.V);
disp('node =');
disp(subst.x);
