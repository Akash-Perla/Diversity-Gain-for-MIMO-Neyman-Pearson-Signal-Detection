% ------------------------------Fig 3-------------------------------------%
% -------- Gaussian reflections and log-normal clutter plus noise --------%


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



%-----------------------------------------------------%
%-----Case: N = 2, M = 2, Q = 2, orth, logn noise-----%
%-----------------------------------------------------%



% Initialize arrays to store results for Q = 2
T_count_array_Q2_orth = zeros(1, length(gamma1_dB)); % Store counts of T > delta_2 for Q = 2
num_simulations = 100; % Number of Monte Carlo simulations
delta_2_Q2_orth = 0.35;

% Loop through each SCNR value for Q = 2
for idx = 1:length(gamma1_linear)
    gamma = gamma1_linear(idx);  % Current SCNR in linear scale

    % Initialize counter for T values greater than delta_2
    T_count = 0;

    for sim = 1:num_simulations
        % Generate uniformly distributed scatterers for Q = 2
        Q = 2; % Number of scatterers
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

        % Generate lognormal noise (w) based on a Gaussian distribution
        mu_logn = -0.45; % Mean of the underlying Gaussian
        sigma_logn = 0.67; % Standard deviation of the underlying Gaussian
        
        % Generate complex lognormal noise
        real_part = exp(mu_logn + sigma_logn * randn(num_rx * num_tx, 1));
        imag_part = exp(mu_logn + sigma_logn * randn(num_rx * num_tx, 1));
        w = (real_part + 1j * imag_part) / sqrt(2); % Normalize noise

        E = 2;
        M = 2;
        scaling_factor = sqrt(E / M);
        C = blkdiag(eye(2), eye(2));
        r = scaling_factor * C * A * zeta + w;

        B = C * A;
        X = gamma * B' * inv(Rw) * B + inv(R);
        
        T = r' * inv(Rw) * B * inv(X) * B' * inv(Rw) * r;

        % Check if T is greater than the threshold
        if T > delta_2_Q2_orth
            T_count = T_count + 1; % Increment counter
        end
    end

    % Store results for this SCNR
    T_count_array_Q2_orth(idx) = T_count; 
    disp(['Q = ', num2str(Q), ' SCNR = ', num2str(gamma1_dB(idx)), ' dB: T count > delta_2 = ', num2str(T_count)]);
end



%--------------------------------------------------------%
%-----Case: N = 2, M = 2, Q = 1500, orth, logn noise-----%
%--------------------------------------------------------%



delta_2_Q1500_orth = 0.35;
% Initialize arrays to store results for Q = 1500
T_count_array_Q1500_orth = zeros(1, length(gamma1_dB)); % Store counts of T > delta_2 for Q = 1500

% Loop through each SCNR value for Q = 1500
for idx = 1:length(gamma1_linear)
    gamma = gamma1_linear(idx);  % Current SCNR in linear scale

    % Initialize counter for T values greater than delta_2
    T_count = 0;

    for sim = 1:num_simulations
        % Generate uniformly distributed scatterers for Q = 1500
        Q = 1500; % Number of scatterers
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
    
        % Generate lognormal noise (w) based on a Gaussian distribution
        mu_logn = -0.45; % Mean of the underlying Gaussian
        sigma_logn = 0.67; % Standard deviation of the underlying Gaussian
        
        % Generate complex lognormal noise
        real_part = exp(mu_logn + sigma_logn * randn(num_rx * num_tx, 1));
        imag_part = exp(mu_logn + sigma_logn * randn(num_rx * num_tx, 1));
        w = (real_part + 1j * imag_part) / sqrt(2); % Normalize noise

        % Calculate r = √(E/M) * C * A * zeta + w
        E = 2;
        M = 2;
        scaling_factor = sqrt(E / M);
        C = blkdiag(eye(2), eye(2));
        r = scaling_factor * C * A * zeta + w;

        % Calculate X = gamma * B^H * Rw^{-1} * B + R^{-1} for current SCNR
        B = C * A;
        X = gamma * B' * inv(Rw) * B + inv(R);

        % Calculate T = r^H * Rw^{-1} * B * X^{-1} * B^H * Rw^{-1} * r
        T = r' * inv(Rw) * B * inv(X) * B' * inv(Rw) * r;

        % Check if T is greater than the threshold
        if T > delta_2_Q1500_orth
            T_count = T_count + 1; % Increment counter
        end
    end

    % Store results for this SCNR
    T_count_array_Q1500_orth(idx) = T_count; % Store count of T > delta_2 for Q = 1500
   disp(['Q = ', num2str(Q), ' SCNR = ', num2str(gamma1_dB(idx)), ' dB: T count > delta_2 = ', num2str(T_count)]);
end



%--------------------------------------------------------%
%-----Case: N = 2, M = 2, Q = 1500, iden, logn noise-----%
%--------------------------------------------------------%



delta_2_Q1500_iden = 0.35;
% Initialize arrays to store results for Q = 1500
T_count_array_Q1500_iden = zeros(1, length(gamma1_dB)); % Store counts of T > delta_2 for Q = 1500

% Loop through each SCNR value for Q = 1500
for idx = 1:length(gamma1_linear)
    gamma = gamma1_linear(idx);  % Current SCNR in linear scale

    % Initialize counter for T values greater than delta_2
    T_count = 0;

    for sim = 1:num_simulations
        % Generate uniformly distributed scatterers for Q = 1500
        Q = 1500; % Number of scatterers
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

        % Generate lognormal noise (w) based on a Gaussian distribution
        mu_logn = -0.45; % Mean of the underlying Gaussian
        sigma_logn = 0.67; % Standard deviation of the underlying Gaussian
        
        % Generate complex lognormal noise
        real_part = exp(mu_logn + sigma_logn * randn(num_rx, 1));
        imag_part = exp(mu_logn + sigma_logn * randn(num_rx, 1));
        w = (real_part + 1j * imag_part) / sqrt(2); % Normalize noise

        % Calculate r = √(E/M) * C * A * zeta + w
        E = 2;
        M = 2;
        scaling_factor = sqrt(E / M);
        Sn = [1 1];
        C = blkdiag(Sn,Sn);
        r = scaling_factor * C * A * zeta + w;

        % Calculate X = gamma * B^H * Rw^{-1} * B + R^{-1} for current SCNR
        B = C * A;
        X = gamma * B' * inv(Rw) * B + inv(R);

        % Calculate T = r^H * Rw^{-1} * B * X^{-1} * B^H * Rw^{-1} * r
        T = r' * inv(Rw) * B * inv(X) * B' * inv(Rw) * r;

        % Check if T is greater than the threshold
        if T > delta_2_Q1500_iden
            T_count = T_count + 1; % Increment counter
        end
    end

    % Store results for this SCNR
    T_count_array_Q1500_iden(idx) = T_count; % Store count of T > delta_2 for Q = 1500
   disp(['Q = ', num2str(Q), ' SCNR = ', num2str(gamma1_dB(idx)), ' dB: T count > delta_2 = ', num2str(T_count)]);
end

disp('P_{M} for each SCNR (Q=2, orth):');
disp(gamma1_dB);
disp(T_count_array_Q2_orth);
disp('P_{M} for each SCNR (Q=1500, orth):');
disp(gamma1_dB);
disp(T_count_array_Q1500_orth);
disp('P_{M} for each SCNR (Q=1500, iden):');
disp(gamma1_dB);
disp(T_count_array_Q1500_iden);

% -------- Plot T_count vs SCNR for both Q values --------
figure;
hold on;
plot(gamma1_dB, T_count_array_Q2_orth/num_simulations, '-o', 'LineWidth', 2, 'DisplayName', 'N = 2, M = 2, Q = 2, orth, logn noise');
plot(gamma1_dB, T_count_array_Q1500_orth/num_simulations, '-s', 'LineWidth', 2, 'DisplayName', 'N = 2, M = 2, Q = 1500, orth, logn noise');
plot(gamma1_dB, T_count_array_Q1500_iden/num_simulations, '-s', 'LineWidth', 2, 'DisplayName', 'N = 2, M = 2, Q = 1500, iden, logn noise');
xlabel('SCNR (dB)');
ylabel('P_{M}')
title('P_{M} vs SCNR');
legend show;
grid on;
xlim([min(gamma1_dB), max(gamma1_dB)]);
hold off;
