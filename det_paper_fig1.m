%-------- Fig 1. System setup --------%

% Coordinates for Transmit Antennas (Tx) and Receive Antennas (Rx)
tx_x = [2000, 6000]; % Tx x-coordinates in meters
tx_y = [-2000, -4000]; % Tx y-coordinates in meters
rx_x = [8000, 4000]; % Rx x-coordinates in meters
rx_y = [2000, 0]; % Rx y-coordinates in meters

% Scatterer distribution range
scatterer_x_min = 300;  % minimum x for scatterers in meters
scatterer_x_max = 1000; % maximum x for scatterers in meters
scatterer_y_min = 9400; % minimum y for scatterers in meters
scatterer_y_max = 10500; % maximum y for scatterers in meters

% Generate uniformly distributed scatterers within the specified region
num_scatterers = 1500; % Number of scatterers
scatterers_x = scatterer_x_min + (scatterer_x_max - scatterer_x_min) * rand(1, num_scatterers);
scatterers_y = scatterer_y_min + (scatterer_y_max - scatterer_y_min) * rand(1, num_scatterers);

% Plot
figure;
hold on;
axis equal;
grid on;

% Plot Transmit Antennas (Tx)
plot(tx_x, tx_y, 'ks', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
text(tx_x(1), tx_y(1), 'Tx1', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
text(tx_x(2), tx_y(2), 'Tx2', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

% Plot Receive Antennas (Rx)
plot(rx_x, rx_y, 'k^', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
text(rx_x(1), rx_y(1), 'Rx1', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right'); %--labeling--%
text(rx_x(2), rx_y(2), 'Rx2', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');     

% Plot scatterers as a densely packed area within the specified range
plot(scatterers_x, scatterers_y, 'k*', 'MarkerSize', 2);
xlim([-6000, 12000]);
ylim([-6000, 12000]);
xlabel('y (m)');
ylabel('x (m)');
title('System setup');
legend({'TX', 'RX', 'Target scatterers'}, 'Location', 'northeast');
hold off;
