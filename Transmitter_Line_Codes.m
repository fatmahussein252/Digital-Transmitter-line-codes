%% Clearing workspace
clear all;
close all;

main();

% Description:
%   This function plots line code signals stored in the 'ensemble' matrix for visualization.
%   It generates subplots, each representing a realization, and plots the corresponding line code signal.
%   The title of each subplot is formed by appending the realization number to the provided string 'str'.
function plot_line_codes(ensemble, str, ylim_start, ylim_end)
    figure;

    for i = 1:4
        subplot(4, 1, i);
        plot(ensemble(i, :));
        title(str, num2str(i));
        ylim([ylim_start ylim_end]);
    end

end

% Description:
%   This function calculates the statistical mean of line code signals across different time steps.
%   It computes the mean value of each time step across all realizations and plots the result.
%   The resulting plot is titled with the provided string 'str'.
function statistical_mean(ensemble, str, num_samples, num_realizations, ylim_start, ylim_end)

    for i = 1:num_samples
        st_mean(1, i) = sum(ensemble(:, i)) / num_realizations;
    end
    
    figure;
    subplot(2, 1, 1);
    plot(st_mean);
    title(str);
    ylim([ylim_start ylim_end]);
end

% Description:
%   Calculate the time mean of line code signals across different realizations
%   It computes the mean value of each realization across time and plots the result.
%   The resulting plot is titled with the provided string 'str'.
function time_mean(ensemble, str, num_samples, num_realizations, ylim_start, ylim_end)

    for i = 1:num_realizations
        t_mean(1, i) = sum(ensemble(i, :)) / num_samples;
    end
    subplot(2, 1, 2);
    plot(t_mean);
    title(str);
    ylim([ylim_start ylim_end]);

end


% Description:
%   This function generates line code signals based on the input data, code type, and amplitude parameters.
%   It supports three types of line codes: polarNRZ, unipolar and polarRZ. For unipolar encoding,
% Output Argument:
%   ensemble: Matrix containing generated line code signals (size: num_realizations x num_samples).
function ensemble = generate_line_code(Data, code_type, A, ensemble, num_realizations)
    
    num_samples_per_bit = 7; % Number of samples per bit
    for i = 1:num_realizations
        if strcmp(code_type, 'unipolar')
            Tx = A * Data(i, :);
        else
            Tx = ((2 * Data(i, :)) - 1) * A;
        end

        Tx2 = repmat(Tx, num_samples_per_bit, 1);
        if strcmp(code_type, 'polarRZ')
            zero_duration_start = 6;
            zero_duration_end = 7;
            Tx2(zero_duration_start:zero_duration_end,:) = 0;
        end
        Tx_out = reshape(Tx2, size(Tx2, 1) * size(Tx2, 2), 1);
        start_index = randi([1 num_samples_per_bit], 1, 1);
        Delayed_Tx = Tx_out(start_index:(size(Tx_out)) - (num_samples_per_bit - start_index) - 1);
        ensemble(i, :) = Delayed_Tx;
    end
end

% Description:
%   Entry point for executing the main functionality of the MATLAB script.
%   This function contains the main code logic that will be executed.

function main()
    % Define the parameters
    A = 4;
    num_realizations = 500;
    num_samples = 700;
    num_bits = 101; % Number of bits in each data instance
    
    % Initialize arrays
    ensemble = zeros(num_realizations, num_samples);
    statistical_mean_array = zeros(1, num_samples);
    time_mean_array = zeros(1, num_realizations);
    Data = randi([0 1], num_realizations, num_bits);
    
    %plots configurations
    ylim_start = -5;
    ylim_end = 5;
    
    %%%% Generation of polar NRZ line code %%%
    ensemble = generate_line_code(Data, "polarNRZ", A, ensemble, num_realizations);
    
    %%%% plot 4 realizations of polar NRZ line code %%%
    plot_line_codes(ensemble, ('polar NRZ line code realisation '), ylim_start, ylim_end);
    
    %%%% calculate statistical and time means of polar NRZ line code %%%
    statistical_mean(ensemble, ("statistical mean - polar NRZ line code"), num_samples, num_realizations, ylim_start, ylim_end);
    time_mean(ensemble, ("time mean - polar NRZ line code"), num_samples, num_realizations, ylim_start, ylim_end);
    
    %%%% Generation of Unipolar line code %%%
    ensemble = generate_line_code(Data, "unipolar", A, ensemble, num_realizations);
    
    %%%% plot 4 realizations of Unipolar line code %%%
    plot_line_codes(ensemble, ('Unipolar line code realisation '), ylim_start, ylim_end);
    
    %%%% calculate statistical and time means of Unipolar line code %%%
    statistical_mean(ensemble, ("statistical mean - Unipolar line code"), num_samples, num_realizations, ylim_start, ylim_end);
    time_mean(ensemble, ("time mean - Unipolar line code"), num_samples, num_realizations, ylim_start, ylim_end);
    
    %%%% Generation of polar RZ line code %%%
    ensemble = generate_line_code(Data, "polarRZ", A, ensemble, num_realizations);
    
    %%%% plot 4 realizations of polar RZ line code %%%
    plot_line_codes(ensemble, ('polar RZ line code realisation '), ylim_start, ylim_end);
    
    %%%% calculate statistical and time means of polar RZ line code %%%
    statistical_mean(ensemble, ("statistical mean - polar RZ line code"), num_samples, num_realizations, ylim_start, ylim_end);
    time_mean(ensemble, ("time mean - polar RZ line code"), num_samples, num_realizations, ylim_start, ylim_end);
    
end
