%% Beam

beam = struct;

% Properties
beam.dim = [1.25; 0.02; 0.08]; % m
beam.vol = prod(beam.dim); % m^3
beam.density = 2690; % kg / m^3
beam.mass = beam.density * beam.vol; % kg

beam.E = 6.45 * 10^10; % N / m^2
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
syms mb mc xc Jc x1 x2 k1 k2 l E Iz;

% Generelized coordinates
syms x;

f = x / l;

q = sym('q', [10 1]);
w = sym('w', size(q));
for i = 1:size(w, 1)
    w(i) = f^(i - 1);
end

M = (mb / l) * int(w * transpose(w), x, -l/2, l/2) + mc * subs(w * transpose(w), x, xc) + Jc * subs(diff(w, x) * transpose(diff(w, x)), x, xc);
K = E * (Iz / l) * int(diff(w, x, 2) * transpose(diff(w, x, 2)), x, -l/2, l/2) + k1 * subs(w * transpose(w), x, x1) + k2 * subs(w * transpose(w), x, x2);

%% Substitution

mb = beam.mass;
mc = cyl.mass;
Jc = cyl.J(3) + cyl.mass * cyl.pos(2)^2;
xc = cyl.pos(1);
x1 = string(1).pos(1);
x2 = string(2).pos(1);
k1 = string(1).stiff;
k2 = string(2).stiff;
l = beam.dim(1);
E = beam.E;
Iz = beam.vol * beam.dim(2)^2 / 12;

subst = struct;
subst.M = double(subs(M));
subst.K = double(subs(K));

%% Eigenvalues and eigenvectors

[subst.V, subst.D2] = eig(subst.K, subst.M);
subst.D = (subst.D2).^(1/2);

for i = 1:size(subst.V, 1)
    subst.D(i, i) = subst.D(i, i) / (2 * pi);
end

subst.f = diag(subst.D);

%% Modes

load('resources/mat/modes2DOF.mat');
modes = load('resources/txt/Project_2018_modes.txt');
P = [-0.592, -0.5:0.1:0.5, 0.592];

y = transpose(subst.V) * subs(w);

dom = -l/2:10^-3:l/2;

y_exp = cell(size(modes, 2), 1);
y_the = cell(size(modes, 2), 1);
y_2dof = cell(size(modes2DOF, 2), 1);

for i = 1:size(modes, 2)
    if i <= size(modes2DOF, 2)
        y_2dof{i} = modes2DOF(1, i) + cm(2) + (dom - cm(1)) * modes2DOF(2, i);
        y_2dof{i} = y_2dof{i} * interp1(P, modes(:, i), 0, 'spline') / y_2dof{i}(dom == 0);
    end
    y_the{i} = double(subs(y(i), x, dom));
    y_the{i} = y_the{i} * interp1(P, modes(:, i), 0, 'spline') / y_the{i}(dom == 0);
end

%% Plot

i = 1;

if i == 1
    plot(dom, zeros(size(dom)), P, modes(:, 1), 'o', dom, y_the{1}, dom, y_2dof{1}, P, modes(:, 2), 'o', dom, y_the{2}, dom, y_2dof{2});
    legend('undeformed', 'experimental n1', 'rayleigh-ritz n1', '2DOF n1', 'experimental n2', 'rayleigh-ritz n2', '2DOF n2');
elseif i == 2
    plot(dom, zeros(size(dom)), P, modes(:, 3), 'o', dom, y_the{3});
    legend('undeformed', 'experimental n3', 'rayleigh-ritz n3');
elseif  i == 3
    plot(dom, zeros(size(dom)), P, modes(:, 4), 'o',  dom, y_the{4}, P, modes(:, 5), 'o', dom, y_the{5});
    legend('undeformed', 'experimental n4', 'rayleigh-ritz n4', 'experimental n5', 'rayleigh-ritz n5');
end

xlabel('x [m]');
ylabel('y [m]');

%% Compare

freq = load('resources/txt/Project_2018_freq.txt');
err = subst.f(1:size(freq, 1)) - freq(:, 1);
rel_err = err ./ freq(:, 1);

%% Display

disp('M =');
disp(subst.M);
disp('K =');
disp(subst.K);
disp('freq =');
disp(subst.f);
