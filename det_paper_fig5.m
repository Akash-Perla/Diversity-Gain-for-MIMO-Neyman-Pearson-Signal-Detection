% -------------Fig 4-------------%
% -------- MIMO and SISO --------%

% Coordinates for Transmit Antennas (Tx) and Receive Antennas (Rx)
tx_x = [2000, 6000];  % Tx x-coordinates in meters
tx_y = [-2000, -4000];  % Tx y-coordinates in meters
rx_x = [8000, 4000];  % Rx x-coordinates in meters
rx_y = [2000, 0];  % Rx y-coordinates in meters

% System parameters
c = 3e8;  % Speed of light in m/s
fc = 1e9; % Carrier frequency in Hz
lambda = c / fc;

% Define gamma1 as an array of SCNR values in dB
gamma1_dB = 0:2:20;
gamma1_linear = 10.^(gamma1_dB / 10);

% Set Q = 1500
Q = 1500;

%--------------------------------------------%
%-----Case: N = 2, M = 2, Q = 1500, MIMO-----%
%--------------------------------------------%

% Initialize arrays to store results for Q = 1500
T_count_array_Q1500_mimo = zeros(1, length(gamma1_dB)); % Store counts of T > delta_2 for Q = 1500
num_simulations = 100000; % Number of Monte Carlo simulations
delta_2_Q1500_mimo = 0.7;

% Loop through each SCNR value for Q = 1500
for idx = 1:length(gamma1_linear)
    gamma = gamma1_linear(idx);  % Current SCNR in linear scale

    % Initialize counter for T values greater than delta_2
    T_count = 0;

    for sim = 1:num_simulations
        % Generate uniformly distributed scatterers for Q = 1500
        scatterers_x = 300 + (1000 - 300) * rand(1, Q);
        scatterers_y = 9400 + (10500 - 9400) * rand(1, Q);

        % Calculate distances and time delays
        num_tx = length(tx_x);
        num_rx = length(rx_x);
        d_tm = zeros(num_rx, num_tx, Q);
        d_rn = zeros(num_rx, num_tx, Q);
        tau_nm = zeros(num_rx, num_tx, Q);

        for n = 1:num_rx
            for m = 1:num_tx
                for q = 1:Q
                    % Distance from Tx to scatterer q
                    d_tm(n, m, q) = sqrt((tx_x(m) - scatterers_x(q))^2 + (tx_y(m) - scatterers_y(q))^2);
                    % Distance from scatterer q to Rx
                    d_rn(n, m, q) = sqrt((rx_x(n) - scatterers_x(q))^2 + (rx_y(n) - scatterers_y(q))^2);
                end
            end
        end

        % Matrix A Generation
        A = zeros(num_rx * num_tx, Q);
        idx_A = 1;
        for n = 1:num_rx
            for m = 1:num_tx
                for q = 1:Q
                    tau_nm(n, m, q) = (d_tm(n, m, q) + d_rn(n, m, q)) / c; % Time delay
                    A(idx_A, q) = exp(-1j * 2 * pi * fc * tau_nm(n, m, q));
                end
                idx_A = idx_A + 1;
            end
        end

        % Covariance matrices for reflection coefficients (R) and noise (Rw)
        R = eye(Q);
        Rw = eye(num_rx);

        % Generate complex Gaussian reflection coefficients (zeta) as zero-mean
        zeta = (randn(Q, 1) + 1j * randn(Q, 1)) / sqrt(2);

        % Generate complex Gaussian noise (w) as zero-mean
        w = (randn(num_rx, 1) + 1j * randn(num_rx, 1)) / sqrt(2);

        E = 2;
        M = 2;
        scaling_factor = sqrt(E / M);
        C = [1, 0, 1, 0; 
             0, 1, 0, 1];
        r = scaling_factor * C * A * zeta + w;

        D = eye(num_rx, num_tx);
        T = r' * (D * D') * r;

        % Check if T is greater than the threshold
        if T > delta_2_Q1500_mimo
            T_count = T_count + 1; % Increment counter
        end
    end

    % Store results for this SCNR
    T_count_array_Q1500_mimo(idx) = T_count; 
    disp(['Q = ', num2str(Q), ' SCNR = ', num2str(gamma1_dB(idx)), ' dB: T count > delta_2 = ', num2str(T_count)]);
end

%--------------------------------------------%
%-----Case: N = 2, M = 2, Q = 1500, SISO-----%
%--------------------------------------------%

% Initialize arrays to store results for Q = 1500
T_count_array_Q1500_siso = zeros(1, length(gamma1_dB)); % Store counts of T > delta_2 for Q = 1500
num_simulations = 100000; % Number of Monte Carlo simulations
delta_2_Q1500_siso = 0.35;

% Loop through each SCNR value for Q = 1500
for idx = 1:length(gamma1_linear)
    gamma = gamma1_linear(idx);  % Current SCNR in linear scale

    % Initialize counter for T values greater than delta_2
    T_count = 0;

    for sim = 1:num_simulations
        % Generate uniformly distributed scatterers for Q = 1500
        scatterers_x = 300 + (1000 - 300) * rand(1, Q);
        scatterers_y = 9400 + (10500 - 9400) * rand(1, Q);

        % Calculate distances and time delays
        num_tx = length(tx_x);
        num_rx = length(rx_x);
        d_tm = zeros(num_rx, num_tx, Q);
        d_rn = zeros(num_rx, num_tx, Q);
        tau_nm = zeros(num_rx, num_tx, Q);

        for n = 1:num_rx
            for m = 1:num_tx
                for q = 1:Q
                    % Distance from Tx to scatterer q
                    d_tm(n, m, q) = sqrt((tx_x(m) - scatterers_x(q))^2 + (tx_y(m) - scatterers_y(q))^2);
                    % Distance from scatterer q to Rx
                    d_rn(n, m, q) = sqrt((rx_x(n) - scatterers_x(q))^2 + (rx_y(n) - scatterers_y(q))^2);
                end
            end
        end

        % Matrix A Generation
        A = zeros(num_rx * num_tx, Q);
        idx_A = 1;
        for n = 1:num_rx
            for m = 1:num_tx
                for q = 1:Q
                    tau_nm(n, m, q) = (d_tm(n, m, q) + d_rn(n, m, q)) / c; % Time delay
                    A(idx_A, q) = exp(-1j * 2 * pi * fc * tau_nm(n, m, q));
                end
                idx_A = idx_A + 1;
            end
        end

        % Covariance matrices for reflection coefficients (R) and noise (Rw)
        R = eye(Q);
        Rw = eye(num_rx * num_tx);

        % Generate complex Gaussian reflection coefficients (zeta) as zero-mean
        zeta = (randn(Q, 1) + 1j * randn(Q, 1)) / sqrt(2);

        % Generate complex Gaussian noise (w) as zero-mean
        w = (randn(num_rx * num_tx, 1) + 1j * randn(num_rx * num_tx, 1)) / sqrt(2);

        E = 2;
        M = 2;
        scaling_factor = sqrt(E / M);
        C = blkdiag(eye(2), eye(2));
        r = scaling_factor * C * A * zeta + w;

        D = [1, 1]; % Row vector representing a^H(theta)
        T = r' * (D * D') * r;

        % Check if T is greater than the threshold
        if T > delta_2_Q1500_siso
            T_count = T_count + 1; % Increment counter
        end
    end

    % Store results for this SCNR
    T_count_array_Q1500_siso(idx) = T_count; 
    disp(['Q = ', num2str(Q), ' SCNR = ', num2str(gamma1_dB(idx)), ' dB: T count > delta_2 = ', num2str(T_count)]);
end

% Plotting the results
figure;
plot(gamma1_dB, T_count_array_Q1500_mimo / num_simulations, '-o', 'DisplayName', 'N = 2, M = 2, Q = 1500, MIMO');
hold on;
plot(gamma1_dB, T_count_array_Q1500_siso / num_simulations, '-x', 'DisplayName', 'N = 2, M = 2, Q = 1500, phased array');
xlabel('SCNR (dB)');
ylabel('Probability (P_M)');
title('Probability P_M vs. SCNR for Q = 1500');
legend('show');
grid on;
