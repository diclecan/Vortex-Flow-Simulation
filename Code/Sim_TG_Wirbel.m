clearvars; close all; clc

%%  PARAMETER
nu = 0.01;
T = 1;
CFLx = 0.2;
u_max = 1;

% Für Netzkonvergenz
N_Werte = [32 64 128 256];

% Für Berichtplots
N_plot = 128;                     % Auflösung für qualitative/zeitliche Plots
nu_values = [0.001 0.01 0.1];     % verschiedene Viskositäten für Vergleich

%%  1) NETZKONVERGENZ
errors = zeros(length(N_Werte),1);
h_Werte = zeros(length(N_Werte),1);

for k = 1:length(N_Werte)

    N = N_Werte(k);

    % Simulation starten
    s = Simulation(N, nu, CFLx, u_max);

    % Gitter
    x = linspace(0, 2*pi, N+1);
    x(end) = [];
    [X, Y] = meshgrid(x);

    X_vektor = X(:);
    Y_vektor = Y(:);

    % Anfangsbedingung Taylor-Green
    omega_start = 2*sin(X_vektor).*sin(Y_vektor);
    omega_vektor = omega_start;

    % Exakte Schrittzahl
    nsteps = ceil(T/s.ht);
    dt = T/nsteps;
    s.ht = dt;

    % Zeitintegration
    for n = 1:nsteps
        Psi = s.PoissonSolver(omega_vektor);
        [u,v] = s.velocities(Psi);
        omega_vektor = s.wirbetransport(u,v,omega_vektor);
    end

    % Analytische Lösung am Endzeitpunkt
    omega_a_vektor = omega_start * exp(-2*nu*T);

    % L2-Fehler
    errors(k) = sqrt(mean((omega_vektor - omega_a_vektor).^2));
    h_Werte(k) = s.hx;

    fprintf("N = %d  Fehler = %.5e\n", N, errors(k));
end

% Gesamtordnung
p = polyfit(log(h_Werte), log(errors), 1);
order = p(1);

fprintf("\nOrdnung = %.4f\n", order);

% Plot Netzkonvergenz
figure
loglog(h_Werte, errors, 'o-', 'LineWidth', 2)
hold on
grid on

ref1 = errors(1) * (h_Werte / h_Werte(1)).^1;
ref2 = errors(1) * (h_Werte / h_Werte(1)).^2;

loglog(h_Werte, ref1, '--', 'LineWidth', 1.5)
loglog(h_Werte, ref2, '-.', 'LineWidth', 1.5)

xlabel('Gitterweite h')
ylabel('L2-Fehler')
title('Netzkonvergenz Taylor-Green-Wirbel')
legend('Numerischer Fehler', 'O(h)', 'O(h^2)', 'Location', 'southwest')
set(gca,'XDir','reverse')

exportgraphics(gcf,'Netzkonvergenz.png','Resolution',300)

%%  2) ANALYTISCHE LÖSUNG:
%     omega(x,y,0) und omega(x,y,T)
N = N_plot;
nu_ana = 0.2;

x = linspace(0, 2*pi, N+1);
x(end) = [];
[X, Y] = meshgrid(x);

X_vektor = X(:);
Y_vektor = Y(:);

omega_start = 2*sin(X_vektor).*sin(Y_vektor);
omega_ana_T = omega_start * exp(-2*nu_ana*T);

figure
contourf(X, Y, reshape(omega_start, N, N), 20, 'LineColor', 'none')
colorbar
clim([-2 2])
axis equal tight
xlabel('x')
ylabel('y')
title('\omega_{ana}(x,y,0), \nu = 0.2')

exportgraphics(gcf,'analytische_Lösung_t0.png','Resolution',300)

figure
contourf(X, Y, reshape(omega_ana_T, N, N), 20, 'LineColor', 'none')
colorbar
clim([-2 2])
axis equal tight
xlabel('x')
ylabel('y')
title('\omega_{ana}(x,y,T), \nu = 0.2')

exportgraphics(gcf,'analytische_Lösung_tT.png','Resolution',300)

%%  3) QUALITATIVER VERGLEICH:
%     numerisch / analytisch / Fehlerfeld
%     für verschiedene Viskositäten
figure

for j = 1:length(nu_values)

    nu_j = nu_values(j);
    s = Simulation(N_plot, nu_j, CFLx, u_max);

    x = linspace(0, 2*pi, N_plot+1);
    x(end) = [];
    [X, Y] = meshgrid(x);

    X_vektor = X(:);
    Y_vektor = Y(:);

    omega_start = 2*sin(X_vektor).*sin(Y_vektor);
    omega_vektor = omega_start;

    nsteps = ceil(T/s.ht);
    dt = T/nsteps;
    s.ht = dt;

    for n = 1:nsteps
        Psi = s.PoissonSolver(omega_vektor);
        [u,v] = s.velocities(Psi);
        omega_vektor = s.wirbetransport(u,v,omega_vektor);
    end

    omega_a_vektor = omega_start * exp(-2*nu_j*T);

    omega_num_mat = reshape(omega_vektor, N_plot, N_plot);
    omega_ana_mat = reshape(omega_a_vektor, N_plot, N_plot);
    err_mat = omega_num_mat - omega_ana_mat;

    subplot(length(nu_values), 3, 3*(j-1)+1)
    contourf(X, Y, omega_num_mat, 20, 'LineColor', 'none')
    colorbar
    clim([-2 2])
    axis equal tight
    xlabel('x')
    ylabel('y')
    title(sprintf('\\omega_{num}, \\nu = %.3g', nu_j))

    subplot(length(nu_values), 3, 3*(j-1)+2)
    contourf(X, Y, omega_ana_mat, 20, 'LineColor', 'none')
    colorbar
    clim([-2 2])
    axis equal tight
    xlabel('x')
    ylabel('y')
    title(sprintf('\\omega_{ana}, \\nu = %.3g', nu_j))

    subplot(length(nu_values), 3, 3*(j-1)+3)
    contourf(X, Y, err_mat, 20, 'LineColor', 'none')
    colorbar
    clim([-0.1 0.1])
    axis equal tight
    xlabel('x')
    ylabel('y')
    title(sprintf('Fehler, \\nu = %.3g', nu_j))
end

sgtitle('Qualitativer Vergleich: numerisch / analytisch / Fehlerfeld')

exportgraphics(gcf,'Qualitativer_Vergleich.png','Resolution',300)

%%  4) ENERGIE UND ENSTROPHIE ÜBER DIE ZEIT
%     für eine ausgewählte Viskosität
nu_eval = 0.01;
s = Simulation(N_plot, nu_eval, CFLx, u_max);

x = linspace(0, 2*pi, N_plot+1);
x(end) = [];
[X, Y] = meshgrid(x);

X_vektor = X(:);
Y_vektor = Y(:);

omega_start = 2*sin(X_vektor).*sin(Y_vektor);
omega_vektor = omega_start;

nsteps = ceil(T/s.ht);
dt = T/nsteps;
s.ht = dt;

t_hist = zeros(nsteps+1,1);

E_num_hist = zeros(nsteps+1,1);
Z_num_hist = zeros(nsteps+1,1);

E_ana_hist = zeros(nsteps+1,1);
Z_ana_hist = zeros(nsteps+1,1);

% Anfangswerte bei t = 0
Psi = s.PoissonSolver(omega_vektor);
[u,v] = s.velocities(Psi);

t_hist(1) = 0;
E_num_hist(1) = 0.5 * sum(u.^2 + v.^2) * s.hx * s.hy;
Z_num_hist(1) = 0.5 * sum(omega_vektor.^2) * s.hx * s.hy;

u_ana =  sin(X_vektor).*cos(Y_vektor) * exp(-2*nu_eval*0);
v_ana = -cos(X_vektor).*sin(Y_vektor) * exp(-2*nu_eval*0);
omega_ana = 2*sin(X_vektor).*sin(Y_vektor) * exp(-2*nu_eval*0);

E_ana_hist(1) = 0.5 * sum(u_ana.^2 + v_ana.^2) * s.hx * s.hy;
Z_ana_hist(1) = 0.5 * sum(omega_ana.^2) * s.hx * s.hy;

for n = 1:nsteps

    Psi = s.PoissonSolver(omega_vektor);
    [u,v] = s.velocities(Psi);
    omega_vektor = s.wirbetransport(u,v,omega_vektor);

    t_now = n * dt;
    t_hist(n+1) = t_now;

    % Numerische Größen
    Psi = s.PoissonSolver(omega_vektor);
    [u,v] = s.velocities(Psi);

    E_num_hist(n+1) = 0.5 * sum(u.^2 + v.^2) * s.hx * s.hy;
    Z_num_hist(n+1) = 0.5 * sum(omega_vektor.^2) * s.hx * s.hy;

    % Analytische Größen
    u_ana =  sin(X_vektor).*cos(Y_vektor) * exp(-2*nu_eval*t_now);
    v_ana = -cos(X_vektor).*sin(Y_vektor) * exp(-2*nu_eval*t_now);
    omega_ana = 2*sin(X_vektor).*sin(Y_vektor) * exp(-2*nu_eval*t_now);

    E_ana_hist(n+1) = 0.5 * sum(u_ana.^2 + v_ana.^2) * s.hx * s.hy;
    Z_ana_hist(n+1) = 0.5 * sum(omega_ana.^2) * s.hx * s.hy;
end

% Plot Energie
figure
plot(t_hist, E_num_hist, 'LineWidth', 2)
hold on
plot(t_hist, E_ana_hist, '--', 'LineWidth', 2)
grid on
xlabel('t')
ylabel('kinetische Energie E(t)')
title(sprintf('Kinetische Energie für \\nu = %.3g', nu_eval))
legend('numerisch', 'analytisch', 'Location', 'northeast')

exportgraphics(gcf,'Energie_num_ana.png','Resolution',300)

% Plot Enstrophie
figure
plot(t_hist, Z_num_hist, 'LineWidth', 2)
hold on
plot(t_hist, Z_ana_hist, '--', 'LineWidth', 2)
grid on
xlabel('t')
ylabel('Enstrophie Z(t)')
title(sprintf('Enstrophie für \\nu = %.3g', nu_eval))
legend('numerisch', 'analytisch', 'Location', 'northeast')

exportgraphics(gcf,'Enstrophie_num_ana.png','Resolution',300)

% Plot relativer Zusammenhang E/Z
figure
plot(t_hist, E_num_hist ./ Z_num_hist, 'LineWidth', 2)
hold on
plot(t_hist, E_ana_hist ./ Z_ana_hist, '--', 'LineWidth', 2)
grid on
xlabel('t')
ylabel('E(t) / Z(t)')
title(sprintf('Verhältnis von Energie und Enstrophie für \\nu = %.3g', nu_eval))
legend('numerisch', 'analytisch', 'Location', 'best')

exportgraphics(gcf,'Energie_Enstrophie.png','Resolution',300)

% Fehler in der Enstrophie
figure
semilogy(t_hist, abs(Z_num_hist - Z_ana_hist), 'LineWidth', 2)
grid on
xlabel('t')
ylabel('|Z_{num} - Z_{ana}|')
title(sprintf('Enstrophiefehler für \\nu = %.3g', nu_eval))

exportgraphics(gcf,'Enstrophiefehler.png','Resolution',300)


%%  5) ANIMATION: analytische Lagrange-Partikel als GIF
%     nu = 0 und nu = 0.2 im selben Fenster
nu_anim = [0 0.2];

% Animationsparameter
T_anim = 5;          % Simulations-/Animationszeit
nFrames = 150;       % Anzahl Frames
dt_anim = T_anim / nFrames;

% Partikelanzahl
nPart_x = 25;
nPart_y = 25;

xp0 = linspace(0.05, 2*pi-0.05, nPart_x);
yp0 = linspace(0.05, 2*pi-0.05, nPart_y);
[XP0, YP0] = meshgrid(xp0, yp0);

% Farben für die Reihen
row_colors = lines(nPart_y);

% Partikel für beide Viskositäten vorbereiten
xp_all = cell(length(nu_anim),1);
yp_all = cell(length(nu_anim),1);
hRows_all = cell(length(nu_anim),1);

for j = 1:length(nu_anim)
    xp_all{j} = XP0(:);
    yp_all{j} = YP0(:);
end

gifname = 'Lagrange_Partikel_vergleich.gif';

fig1 = figure('Name', 'Analytische Lagrange-Partikel: Vergleich');
clf

for j = 1:length(nu_anim)

    subplot(1,2,j)
    hold on

    hRows = gobjects(nPart_y,1);
    xp = xp_all{j};
    yp = yp_all{j};

    for r = 1:nPart_y
        idx = (r-1)*nPart_x + (1:nPart_x);
        hRows(r) = plot(xp(idx), yp(idx), '.', ...
            'Color', row_colors(r,:), 'MarkerSize', 10);
    end

    axis equal
    axis([0 2*pi 0 2*pi])
    xlabel('x')
    ylabel('y')
    title(sprintf('\\nu = %.3g, t = %.3f', nu_anim(j), 0))

    hRows_all{j} = hRows;
end

sgtitle('Analytische Lagrange-Partikel')
drawnow

frame = getframe(fig1);
im = frame2im(frame);
[imind, cm] = rgb2ind(im, 256);
imwrite(imind, cm, gifname, 'gif', 'Loopcount', inf, 'DelayTime', T_anim/nFrames);

for n = 1:nFrames

    t_now = (n-1) * dt_anim;
    t_mid = t_now + 0.5*dt_anim;

    for j = 1:length(nu_anim)

        nu_j = nu_anim(j);
        xp = xp_all{j};
        yp = yp_all{j};

        u1 =  sin(xp).*cos(yp) * exp(-2*nu_j*t_now);
        v1 = -cos(xp).*sin(yp) * exp(-2*nu_j*t_now);

        xm = xp + 0.5*dt_anim*u1;
        ym = yp + 0.5*dt_anim*v1;

        u2 =  sin(xm).*cos(ym) * exp(-2*nu_j*t_mid);
        v2 = -cos(xm).*sin(ym) * exp(-2*nu_j*t_mid);

        xp = xp + dt_anim*u2;
        yp = yp + dt_anim*v2;

        % periodische Randbedingungen
        xp = mod(xp, 2*pi);
        yp = mod(yp, 2*pi);

        xp_all{j} = xp;
        yp_all{j} = yp;

        % Plotdaten updaten
        hRows = hRows_all{j};
        for r = 1:nPart_y
            idx = (r-1)*nPart_x + (1:nPart_x);
            set(hRows(r), 'XData', xp(idx), 'YData', yp(idx));
        end

        subplot(1,2,j)
        title(sprintf('\\nu = %.3g, t = %.3f', nu_j, n*dt_anim))
    end

    drawnow limitrate

    frame = getframe(fig1);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    imwrite(imind, cm, gifname, 'gif', 'WriteMode', 'append', ...
        'DelayTime', T_anim/nFrames);
end

%%  6) ANIMATION: analytische 3D-Stromfunktion als GIF
%     nu = 0 und nu = 0.2 im selben Fenster
nu_anim = [0 0.2];

% Gitter für die analytische Stromfunktion
N_anim = 64;
x = linspace(0, 2*pi, N_anim+1);
x(end) = [];
[X, Y] = meshgrid(x);

% Animationsparameter
T_anim = 5;
nFrames = 150;
dt_anim = T_anim / nFrames;

gifname = 'Stromfunktion_3D_vergleich.gif';

fig2 = figure('Name', 'Analytische 3D-Stromfunktion: Vergleich');
clf

hSurf_all = gobjects(length(nu_anim),1);
hCont_all = gobjects(length(nu_anim),1);

for j = 1:length(nu_anim)

    nu_j = nu_anim(j);
    Psi = sin(X).*sin(Y) * exp(-2*nu_j*0);

    subplot(1,2,j)
    hSurf_all(j) = surf(X, Y, Psi, 'EdgeColor', 'none');
    hold on
    [~, hCont_all(j)] = contour3(X, Y, Psi, 18, 'k');

    xlabel('x')
    ylabel('y')
    zlabel('\psi')
    title(sprintf('\\nu = %.3g, t = %.3f', nu_j, 0))
    view(45,30)
    axis tight
    zlim([-1 1])
end

sgtitle('Analytische 3D-Stromfunktion')
drawnow

frame = getframe(fig2);
im = frame2im(frame);
[imind, cm] = rgb2ind(im, 256);
imwrite(imind, cm, gifname, 'gif', 'Loopcount', inf, 'DelayTime', T_anim/nFrames);

for n = 1:nFrames

    t_now = n * dt_anim;

    for j = 1:length(nu_anim)

        nu_j = nu_anim(j);
        Psi = sin(X).*sin(Y) * exp(-2*nu_j*t_now);

        subplot(1,2,j)

        set(hSurf_all(j), 'ZData', Psi)

        delete(hCont_all(j))
        [~, hCont_all(j)] = contour3(X, Y, Psi, 18, 'k');

        title(sprintf('\\nu = %.3g, t = %.3f', nu_j, t_now))
        view(45,30)
        zlim([-1 1])
    end

    drawnow limitrate

    frame = getframe(fig2);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    imwrite(imind, cm, gifname, 'gif', 'WriteMode', 'append', ...
        'DelayTime', T_anim/nFrames);
end