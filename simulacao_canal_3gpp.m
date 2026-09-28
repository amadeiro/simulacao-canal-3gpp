clear; clc; close all;

rng(42);

forcar_LOS = [];

fc = 3e9;
c = 3e8;
wavelength = c / fc;
N = 100;
vrx = 5.0;

BS_pos = [0, 0, 25];
UT_pos = [100, 100, 1.5];

d2D = norm(UT_pos(1:2) - BS_pos(1:2));
d3D = norm(UT_pos - BS_pos);
v_hat = [1, 0, 0];
v_vec = vrx * v_hat;

phi_LoS = atan2d(BS_pos(2) - UT_pos(2), BS_pos(1) - UT_pos(1));
elev_LoS = atan2d(BS_pos(3) - UT_pos(3), d2D);
theta_LoS = 90 - elev_LoS;

h_UT = UT_pos(3);
if h_UT <= 13
    C_hUT = 0;
else
    C_hUT = ((h_UT - 13) / 10)^1.5;
end

if d2D <= 18.0
    P_LoS = 1.0;
else
    P_LoS = ((18.0 / d2D) + exp(-d2D / 63.0) * (1.0 - 18.0 / d2D)) * ...
            (1 + C_hUT * (5/4) * (d2D / 100)^3 * exp(-d2D / 150));
end

if isempty(forcar_LOS)
    is_LOS = (rand() < P_LoS);
else
    is_LOS = logical(forcar_LOS);
end

log_fc = log10(max(fc / 1e9, 6));

if is_LOS
    mu_lgDS = -6.955 - 0.0963 * log_fc;
    sigma_lgDS = 0.66;

    mu_lgASA = 1.81;
    sigma_lgASA = 0.20;

    mu_lgZSA = 0.95;
    sigma_lgZSA = 0.16;

    mu_K = 9.0; sigma_K = 3.5;
    r_tau = 2.5;
else
    mu_lgDS = -6.28 - 0.204 * log_fc;
    sigma_lgDS = 0.39;

    mu_lgASA = 2.08 - 0.27 * log_fc;
    sigma_lgASA = 0.11;

    mu_lgZSA = 1.512 - 0.3236 * log_fc;
    sigma_lgZSA = 0.16;

    mu_K = 0.0; sigma_K = 0.0;
    r_tau = 2.3;
end
zeta = 3.0;

sigma_tau = 10^(mu_lgDS + sigma_lgDS * randn());
sigma_ASA = min(10^(mu_lgASA + sigma_lgASA * randn()), 104.0);
sigma_ZSA = min(10^(mu_lgZSA + sigma_lgZSA * randn()), 52.0);

if is_LOS
    K_dB = mu_K + sigma_K * randn();
    K_lin = 10^(K_dB / 10.0);
else
    K_lin = 0.0;
end

mu_tau = r_tau * sigma_tau;
tau_raw = -mu_tau * log(rand(1, N));
tau_prime = tau_raw - min(tau_raw);
tau = sort(tau_prime);

xi = zeta * randn(1, N);
P_hat = exp(-tau .* (r_tau - 1.0) ./ (r_tau * sigma_tau)) .* 10.^(-xi ./ 10.0);

P = zeros(1, N);
if is_LOS
    P(1) = K_lin / (K_lin + 1.0);
    Omega_c = sum(P_hat(2:end));
    P(2:end) = (1.0 / (K_lin + 1.0)) .* (P_hat(2:end) ./ Omega_c);
else
    Omega_c = sum(P_hat);
    P = P_hat ./ Omega_c;
end

alpha = sqrt(P);

tau_medio = sum(P .* tau);
sigma_tau_obtido = sqrt(sum(P .* tau.^2) - tau_medio^2);

p_ratio = max(P ./ max(P), 1e-12);

phi_double_prime = 1.42 * sigma_ASA * sqrt(-log(p_ratio));
U = randi([0, 1], 1, N)*2 - 1;
Y = (sigma_ASA / 7.0) * randn(1, N);
phi_prime = U .* phi_double_prime + Y + phi_LoS;
if is_LOS
    phi_prime(1) = phi_LoS;
end
phi_n = mod(phi_prime, 360.0);

theta_double_prime = -sigma_ZSA * log(p_ratio);
U_z = randi([0, 1], 1, N)*2 - 1;
Y_z = (sigma_ZSA / 7.0) * randn(1, N);
theta_prime = U_z .* theta_double_prime + Y_z + theta_LoS;
if is_LOS
    theta_prime(1) = theta_LoS;
end
theta_n = mod(theta_prime, 360.0);
theta_n(theta_n > 180) = 360 - theta_n(theta_n > 180);

phi_rad = deg2rad(phi_n);
theta_rad = deg2rad(theta_n);

r_vecs = [cos(phi_rad).*sin(theta_rad); ...
          sin(phi_rad).*sin(theta_rad); ...
          cos(theta_rad)];

nu = (vrx / wavelength) * (v_hat * r_vecs);
phi_bar = 2 * pi * (fc + nu) .* tau;

rho_lim = 0.95;

kappa_vec = [0, logspace(2, 7, 5000)];
rho_freq = abs(P * exp(-1j * 2 * pi * tau.' * kappa_vec));
Bc = limiar_correlacao(kappa_vec, rho_freq, rho_lim);

sigma_time_vec = [0, logspace(-6, log10(0.05), 20000)];
v_list = [5, 50];
rho_time = zeros(numel(v_list), numel(sigma_time_vec));
Tc = zeros(1, numel(v_list));
for iv = 1:numel(v_list)
    nu_v = nu * (v_list(iv) / vrx);
    rho_time(iv, :) = abs(P * exp(1j * 2 * pi * nu_v.' * sigma_time_vec));
    Tc(iv) = limiar_correlacao(sigma_time_vec, rho_time(iv, :), rho_lim);
end

nu_medio = sum(P .* nu);
sigma_nu = sqrt(sum(P .* nu.^2) - nu_medio^2);
Tc_aprox = @(s_nu) sqrt(2 * (1 - rho_lim)) ./ (2 * pi * s_nu);

v_sweep = logspace(0, 2, 40);
fD_sweep = v_sweep / wavelength;
Tc_sweep = zeros(size(v_sweep));
for iv = 1:numel(v_sweep)
    nu_v = nu * (v_sweep(iv) / vrx);
    Tc_sweep(iv) = limiar_correlacao(sigma_time_vec, ...
        abs(P * exp(1j * 2 * pi * nu_v.' * sigma_time_vec)), rho_lim);
end
sigma_nu_sweep = sigma_nu * (v_sweep / vrx);

fprintf('===============================================================\n');
fprintf('  RESULTADOS DA SIMULAÇÃO DE CANAL (3GPP TR 38.901 UMa - MATLAB)\n');
fprintf('===============================================================\n');
fprintf(' Frequência portadora : %.1f GHz\n', fc / 1e9);
fprintf(' Probabilidade LoS    : %.3f\n', P_LoS);
if is_LOS
    fprintf(' Condição de Enlace   : LOS (Com Visada Direta)\n');
else
    fprintf(' Condição de Enlace   : NLOS (Sem Visada Direta)\n');
end
fprintf(' Distância 2D / 3D    : %.2f m / %.2f m\n', d2D, d3D);
fprintf(' Ângulos LoS (UT)     : phi = %.2f° | elev = %.2f° | zenital = %.2f°\n', phi_LoS, elev_LoS, theta_LoS);
fprintf(' Delay Spread sorteado: %.2f ns\n', sigma_tau * 1e9);
fprintf(' Delay Spread obtido  : %.2f ns (RMS do PDP gerado)\n', sigma_tau_obtido * 1e9);
fprintf(' Espalhamento Azimute : %.2f°\n', sigma_ASA);
fprintf(' Espalhamento Zenital : %.2f°\n', sigma_ZSA);
if is_LOS
    fprintf(' Fator K de Rice      : %.2f dB (%.2f linear)\n', K_dB, K_lin);
else
    fprintf(' Fator K de Rice      : 0 (NLOS)\n');
end
fprintf(' Número de Componentes: %d\n', N);
fprintf(' Doppler máximo       : %.2f Hz (v = %.1f m/s)\n', vrx / wavelength, vrx);
fprintf(' Doppler médio / RMS  : %.2f Hz / %.2f Hz (v = %.1f m/s)\n', nu_medio, sigma_nu, vrx);
fprintf(' Banda de Coerência   : %s (|rho| >= %.2f)\n', formatar(Bc, 1e3, 'kHz', kappa_vec(end)), rho_lim);
for iv = 1:numel(v_list)
    fprintf(' Tempo de Coerência   : %s (v = %d m/s) | aprox. 0.0503/sigma_nu = %.2f ms\n', ...
        formatar(Tc(iv), 1e-3, 'ms', sigma_time_vec(end)), v_list(iv), ...
        Tc_aprox(sigma_nu * v_list(iv) / vrx) * 1e3);
end
fprintf('===============================================================\n\n');

nova_figura('Perfil de Atraso e Potencia (PDP)');
stem(tau * 1e6, P, 'filled', 'LineWidth', 1.2, 'MarkerSize', 4);
set(gca, 'YScale', 'log');
grid on;
xlabel('Domínio de Atraso \tau (\mus)', 'FontSize', 11);
ylabel('PDP  \alpha_n^2', 'FontSize', 11);
title(sprintf('Figura 1: Perfil de Atraso e Potência (\\sigma_\\tau = %.2f ns)', sigma_tau * 1e9), 'FontSize', 12);

nova_figura('PAS Azimute');
plot_pas(phi_n, P);
title(sprintf('Figura 2: Espectro Angular de Potência — Azimute (\\sigma_{\\phi,AoA} = %.2f°)', sigma_ASA), 'FontSize', 12);

nova_figura('PAS Elevacao');
plot_pas(theta_n, P);
title(sprintf('Figura 3: Espectro Angular de Potência — Zenital (\\sigma_{\\theta,AoA} = %.2f°)', sigma_ZSA), 'FontSize', 12);

nova_figura('Vetores Direcao de Chegada');
P_dB = 10 * log10(P);
quiver3(zeros(1, N), zeros(1, N), zeros(1, N), r_vecs(1, :), r_vecs(2, :), r_vecs(3, :), 0, ...
        'Color', [0.6 0.6 0.6], 'LineWidth', 0.6, 'ShowArrowHead', 'off'); hold on;
scatter3(r_vecs(1, :), r_vecs(2, :), r_vecs(3, :), 30, P_dB, 'filled');
if is_LOS
    quiver3(0, 0, 0, r_vecs(1, 1), r_vecs(2, 1), r_vecs(3, 1), 0, 'r', 'LineWidth', 2.5, 'MaxHeadSize', 0.3);
end
quiver3(0, 0, 0, v_hat(1), v_hat(2), v_hat(3), 0, 'k', 'LineWidth', 2, 'MaxHeadSize', 0.3);
cb = colorbar; cb.Label.String = 'Potência \alpha_n^2 (dB)';
axis equal; grid on; view(35, 25);
xlim([-1 1]); ylim([-1 1]); zlim([-1 1]);
xlabel('x'); ylabel('y'); zlabel('z');
if is_LOS
    legend('Vetores r_n', 'Componentes (cor = potência)', 'Componente LoS', 'Direção v_{rx}', 'Location', 'bestoutside');
else
    legend('Vetores r_n', 'Componentes (cor = potência)', 'Direção v_{rx}', 'Location', 'bestoutside');
end
title('Figura 4: Vetores Direção das Componentes Multipercurso', 'FontSize', 12);

nova_figura('Espectro Doppler');
stem(nu, P, 'filled', 'Color', [0.85 0.33 0.1], 'LineWidth', 1.2, 'MarkerSize', 4);
set(gca, 'YScale', 'log');
grid on;
xlim(1.1 * [-1 1] * vrx / wavelength);
xlabel('Domínio de Desvio Doppler \nu (Hz)', 'FontSize', 11);
ylabel('Espectro Doppler  \alpha_n^2', 'FontSize', 11);
title('Figura 5: Espectro Doppler das Componentes Multipercurso', 'FontSize', 12);

delta_t_ex = 1e-7;
Nt = 1e5;
t_ex = linspace(0, 5 * delta_t_ex, Nt);
s_ex = gerar_pulso(t_ex, 0, delta_t_ex);
f_ex = linspace(-5 / delta_t_ex, 5 / delta_t_ex, 4001);
S_ex = abs(delta_t_ex * sinc_local(f_ex * delta_t_ex));

nova_figura('Pulso Transmitido');
subplot(2, 1, 1);
plot(t_ex * 1e6, s_ex, 'k', 'LineWidth', 1.5); grid on;
ylim([-0.1 1.2]);
xlabel('Tempo (\mus)'); ylabel('s(t)');
title(sprintf('Figura 6: Pulso Transmitido (\\delta t = %g ns)', delta_t_ex * 1e9));
subplot(2, 1, 2);
plot(f_ex / 1e6, S_ex / max(S_ex), 'b', 'LineWidth', 1.5); hold on; grid on;
xline([-1 1] / delta_t_ex / 1e6, 'r--', 'B_w \approx 1/\delta t', 'LineWidth', 1.2);
xlabel('Frequência (MHz)'); ylabel('|S(f)| normalizado');
title('Espectro de Magnitude do Pulso');

delta_t_list = [1e-7, 1e-5, 1e-4];
for k = 1:numel(delta_t_list)
    dt = delta_t_list(k);
    t_max = max(5 * dt, max(tau) + 2 * dt);
    t_vec = linspace(0, t_max, Nt);
    s_tx = gerar_pulso(t_vec, 0, dt);
    r_t = sinal_recebido(t_vec, tau, alpha, nu, phi_bar, dt);

    nova_figura(sprintf('Sinal Recebido dt = %g s', dt));
    plot(t_vec * 1e6, s_tx, 'k--', 'LineWidth', 1.5); hold on;
    plot(t_vec * 1e6, abs(r_t), 'r', 'LineWidth', 1.2);
    grid on; legend('Pulso Transmitido |s(t)|', 'Sinal Recebido |r(t)|');
    xlabel('Tempo (\mus)'); ylabel('Amplitude');
    title(sprintf('Figura %d: Sinal recebido, \\delta t = %g s, \\sigma_\\tau = %.2f ns', ...
        6 + k, dt, sigma_tau * 1e9), 'FontSize', 12);
end

nova_figura('Banda de Coerencia Bc');
semilogx(kappa_vec(2:end) / 1e3, rho_freq(2:end), 'LineWidth', 2, 'Color', [0 0.5 0]); hold on;
yline(rho_lim, 'r--', sprintf('\\rho_B = %.2f', rho_lim), 'LineWidth', 1.5);
if ~isnan(Bc)
    xline(Bc / 1e3, 'k:', sprintf('B_c = %.1f kHz', Bc / 1e3), 'LineWidth', 1.5, ...
          'LabelVerticalAlignment', 'bottom');
end
grid on; ylim([0 1.05]);
xlabel('Desvio de Frequência \kappa (kHz) — escala log', 'FontSize', 11);
ylabel('|\rho_{TT}(\kappa; 0)|', 'FontSize', 11);
title('Figura 10: Autocorrelação em Frequência e Banda de Coerência B_c', 'FontSize', 12);

nova_figura('Tempo de Coerencia Tc');
estilos = {'b', 'r--'};
leg = cell(1, numel(v_list));
for iv = 1:numel(v_list)
    semilogx(sigma_time_vec(2:end) * 1e3, rho_time(iv, 2:end), estilos{iv}, 'LineWidth', 2); hold on;
    leg{iv} = sprintf('v_{rx} = %d m/s (%d km/h), T_c = %s', v_list(iv), round(v_list(iv) * 3.6), ...
                      formatar(Tc(iv), 1e-3, 'ms', sigma_time_vec(end)));
end
yline(rho_lim, 'k:', sprintf('\\rho_T = %.2f', rho_lim), 'LineWidth', 1.2, 'HandleVisibility', 'off');
for iv = 1:numel(v_list)
    if ~isnan(Tc(iv))
        xline(Tc(iv) * 1e3, estilos{iv}(1), 'LineWidth', 1, 'LineStyle', ':', 'HandleVisibility', 'off');
    end
end
grid on; ylim([0 1.05]);
legend(leg, 'Location', 'southwest');
xlabel('Desvio de Tempo \sigma (ms) — escala log', 'FontSize', 11);
ylabel('|\rho_{TT}(0; \sigma)|', 'FontSize', 11);
title('Figura 11: Autocorrelação Temporal e Tempo de Coerência T_c vs Velocidade', 'FontSize', 12);

nova_figura('Tc vs Doppler');
loglog(fD_sweep, Tc_sweep * 1e3, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 5); hold on;
loglog(fD_sweep, Tc_aprox(sigma_nu_sweep) * 1e3, 'r-', 'LineWidth', 1.8);
for iv = 1:numel(v_list)
    loglog(v_list(iv) / wavelength, Tc(iv) * 1e3, 'ks', 'MarkerSize', 12, 'LineWidth', 1.8, ...
           'HandleVisibility', 'off');
    text(v_list(iv) / wavelength * 1.12, Tc(iv) * 1e3 * 1.12, sprintf('%d m/s', v_list(iv)), 'FontSize', 10);
end
grid on;
legend('T_c simulado (|\rho_{TT}(0;\sigma)| = 0.95)', ...
       sprintf('Aproximação T_c \\approx 0.0503/\\sigma_\\nu  (\\sigma_\\nu = %.2f f_D)', sigma_nu / (vrx / wavelength)), ...
       'Location', 'northeast');
xlabel('Doppler máximo f_D = v_{rx}/\lambda (Hz)', 'FontSize', 11);
ylabel('Tempo de Coerência T_c (ms)', 'FontSize', 11);
title('Figura 12: Tempo de Coerência vs Desvio Doppler (T_c \propto 1/f_D)', 'FontSize', 12);

function fig = nova_figura(nome)
    fig = figure('Name', nome, 'Color', 'w');
    if exist('theme', 'file')
        theme(fig, 'light');
    end
end

function s = gerar_pulso(t, t0, delta_t)
    s = double(t >= t0 & t < t0 + delta_t);
end

function r_t = sinal_recebido(t_vec, tau, alpha, nu, phi_bar, delta_t)
    r_t = zeros(size(t_vec));
    for n = 1:length(tau)
        fase_t = phi_bar(n) - 2 * pi * nu(n) * t_vec;
        r_t = r_t + alpha(n) * exp(-1j * fase_t) .* gerar_pulso(t_vec, tau(n), delta_t);
    end
end

function x_c = limiar_correlacao(x, rho, limiar)
    idx = find(rho < limiar, 1, 'first');
    if isempty(idx)
        x_c = NaN;
    else
        x_c = x(idx - 1) + (limiar - rho(idx - 1)) * (x(idx) - x(idx - 1)) / (rho(idx) - rho(idx - 1));
    end
end

function txt = formatar(valor, escala, unidade, x_max)
    if isnan(valor)
        txt = sprintf('> %.1f %s', x_max / escala, unidade);
    else
        txt = sprintf('%.2f %s', valor / escala, unidade);
    end
end

function y = sinc_local(x)
    y = ones(size(x));
    nz = x ~= 0;
    y(nz) = sin(pi * x(nz)) ./ (pi * x(nz));
end

function plot_pas(ang_deg, P)
    P_dB = 10 * log10(P);
    piso = 10 * floor(min(P_dB) / 10);
    polarplot(deg2rad(ang_deg), P_dB - piso, 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 4);
    pax = gca;
    ticks = 0:10:(-piso);
    pax.RLim = [0, -piso];
    pax.RTick = ticks;
    pax.RTickLabel = arrayfun(@(v) sprintf('%d dB', v + piso), ticks, 'UniformOutput', false);
end
