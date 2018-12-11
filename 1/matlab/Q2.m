Q1;
close;

%% Parameters

g = 9.81;
m = 4.484;
k = 10500;

%%%%% Computing

d = 100;
r = 3;

%% 5

%%%%% Load

FRF = load('resources/txt/FRF_1DOF.txt');

%%%%% Bode

X = (FRF(:, 3).^2 + FRF(:, 2).^2).^(1/2);
phi = atan2(FRF(:, 3), FRF(:, 2));

subplot(2,1,1);
plot(FRF(:, 1) * 2 * pi, X);
subplot(2,1,2);
plot(FRF(:, 1) * 2 * pi, phi);

%%%%% Interpolation around the peak

[~, i] = max(X);

interval = linspace(FRF(i - r, 1), FRF(i + r, 1), d);

interp = interp1(FRF(i - r:i + r, 1), X(i - r:i + r), interval, 'spline');

[X_max, j] = max(interp);

f = interval(j);

%%%%% Half-power

f_l = interp1(X(1:i - 1), FRF(1:i - 1, 1), X_max / 2^(1/2), 'spline');
f_r = interp1(X(i + 1:end), FRF(i + 1:end, 1), X_max / 2^(1/2), 'spline');

w_bode = f * 2 * pi;

delta_w = (f_r - f_l) * 2 * pi;

%%%%% Quality factor and damping ratio

Q_bode = w_0 / delta_w;

eps_bode = 1 / (2 * Q_bode);

%% 6

%%%%% Nyquist

%plot(FRF(:, 2), FRF(:, 3));
%daspect([1 1 1]);

%%%%% Damping ratio

for i = 1:size(FRF, 1)
    if FRF(i, 2) * FRF(1, 2) < 0
        break;
    end
end

f_res = interp1(FRF(i - 2:i + 1, 2), FRF(i - 2:i + 1, 1), 0, 'spline');
i_res = interp1(FRF(i - 2:i + 1, 1), FRF(i - 2:i + 1,3), f_res, 'spline');

w_nyq = f_res * 2 * pi;

Q_nyq = abs(i_res) * g * m;

eps_nyq = 1 / (2 * Q_nyq);

clearvars -except w_0 w_d eps X_max w_bode Q_bode eps_bode w_nyq Q_nyq eps_nyq;
