%% Parameters

g = 9.81;
m = 4.484;
k = 10500;

%%%%% Computing

s = 15;
d = 100;
n = 50;
r = 5;

%% 1

w_0 = (k / m)^(1/2);

%% 3

%%%%% Load

IRF = load('resources/txt/IRF_1DOF.txt');
IRF = IRF(s:end - s, :);

%%%%% Time response

plot(IRF(:, 1), IRF(:, 2));

%%%%% Compute 1^st and n^th maxima

for i = 2:size(IRF, 1) - 1
    if IRF(i, 2) > IRF(i - 1, 2) && IRF(i, 2) > IRF(i + 1, 2)
        i_1 = i;
        break;
    end
end

l = 1;
i_n = i_1 + 1;

while l < n
    if IRF(i_n, 2) > IRF(i_n - 1, 2) && IRF(i_n, 2) > IRF(i_n + 1, 2)
        l = l + 1;
    end
    i_n = i_n + 1;
end

i_n = i_n - 1;

%%%%% Interpolation around both peaks

interval_1 = linspace(IRF(i_1 - 1, 1), IRF(i_1 + r, 1), d);
interval_n = linspace(IRF(i_n - r, 1), IRF(i_n + 1, 1), d);

interp_1 = interp1(IRF(i_1 - 1:i_1 + r, 1), IRF(i_1 - 1:i_1 + r, 2), interval_1, 'spline');
interp_n = interp1(IRF(i_n - r:i_n + 1, 1), IRF(i_n - r:i_n + 1, 2), interval_n, 'spline');

[a_1, j_1] = max(interp_1);
[a_n, j_n] = max(interp_n);

t_1 = interval_1(j_1);
t_n = interval_n(j_n);

%%%%% Natural frequency and damping ratio

w_d = 2 * pi * (n - 1) / (t_n - t_1);

delta = log(a_1 / a_n) / (n - 1);

eps = delta / (2 * pi);

clearvars -except w_0 w_d eps;
