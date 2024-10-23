% wrote by Ernesto Neira, 
% Directed Energy Research Centre, 
% Technology Innovation Institute, Abu Dhabi, UAE.
% Oct 2024

function [varargout] = matrixPencilMethod(t, signal, varargin)
    % matrixPencilMethod applies the matrix pencil method to extract poles from
    % a signal, then reconstructs the signal and its integral.
    %
    % Syntax:
    %   [poles, residues, x_r, int_x] = matrixPencilMethod(t, signal, 'Accuracy', 0.0001, 'method', 'CNRU')
    %
    % Inputs:
    %   t       - Time vector
    %   signal  - Input signal
    %   varargin- Name-value pairs for optional parameters:
    %             'Accuracy' - Accuracy or tolerance for pole estimation (default: 0.05)
    %             'Method'   - Method for pole filtering ('CNRB' or 'CNRU',
    %             default: 'CNRB'), Bilateral or unilateral.
    %             'DC' - Eliminate the DC of the Signal y or n, (default y)
    %
    % Outputs:
    %   poles   - Extracted poles of the signal
    %   residues- Residues associated with each pole
    %   x_r     - Reconstructed signal
    %   int_x   - Integral of the reconstructed signal
    %
    % Example of usage:
    %   clear all; close all;
    %   t = linspace(0, 1, 1000);
    %   alpha = [-0.9, -0.7, -0.8, -0.1];
    %   f = [1, 15, 25, 20];
    %   c = [10, 2, 1, 5];
    %   x = zeros(1, length(t));
    %   for i = 1:length(f)
    %       x = c(i) * exp(2 * pi * alpha(i) * t) .* cos(2 * pi * f(i) * t) + x;
    %   end
    %   [poles, residues, x_r, int_x] = matrixPencilMethod(t, x, 'Accuracy', 0.0001, 'method', 'CNRU');
    %   plot(t, x, 'k', 'LineWidth', 3); hold on;
    %   plot(t, x_r, 'r'); legend('Original Signal', 'Reconstructed Signal');
    
    % Default options
    options.Method = 'CNRB';  % Default method for filtering poles
    options.Accuracy = 0.05;  % Default accuracy for pole estimation
    options.DC ='n'; %Eliminate the DC of the signal

    % Check that varargin contains name-value pairs
    if mod(length(varargin), 2) ~= 0
        error('The parameter must be pairs');
    end

    % Parse optional parameters to overwrite defaults
    for i = 1:2:length(varargin)
        if isempty(varargin{i+1})
            % If value is empty, use default value
            continue;
        end
        switch lower(varargin{i})
            case 'method'
                options.Method = varargin{i+1};  % Override method
            case 'accuracy'
                options.Accuracy = varargin{i+1};  % Override accuracy
            case 'dc'
                options.DC= varargin{i+1};
            otherwise
                error(['Unknown option: ', varargin{i}]);
        end
    end

    % Remove the mean from the signal
    if options.DC=='n';
        y = signal - mean(signal);
    else
        y=signal;
    end

    % Call matrixPencil to compute the poles of the signal
    poles = matrixPencil(t, y, options.Accuracy, options.Method);

    % Call reconstructMPMint2 to reconstruct the signal and compute the integral
    [residues, x_r, int_x] = reconstructMPMint2(poles, t, signal, options.Method);

    % Return the poles, residues, reconstructed signal, and integral
    variables = {poles, residues, x_r, int_x};
    num_outputs = nargout;  % Number of requested outputs
    
    % Assign the output variables
    for i = 1:num_outputs
        varargout{i} = variables{i};
    end
end

function poles = matrixPencil(t, y, numPoles, Method)
    % matrixPencil applies the matrix pencil method to estimate the poles
    % of the input signal y based on the time vector t and given method.

    % Normalize the signal
    y = y / max(y);
    N = length(y);  % Length of the signal
    K = round(N / 2);  % Parameter K for matrix pencil method

    % Ensure the signal is a column vector
    y = y(:);

    % Create the data matrix A for the matrix pencil method
    A = zeros(N-K, K+1);
    for i = 1:K+1
        A(:, i) = y(i:i+N-K-1);
    end

    % Perform singular value decomposition (SVD)
    [u, s, v] = svd(A(:, 2:K+1), 0);

    % Determine the number of poles (L) based on numPoles input
    if numPoles == 0
        % Auto-rank selection with tolerance
        M = 10;
        if Method == 'CNRU'; M = 20; end
        tol = s(1, 1) / M;
        L = rank(s, tol);
    elseif numPoles > 0 && numPoles < 1
        % If numPoles is a decimal, use it as a precision factor
        if Method == 'CNRU'; numPoles = numPoles / 2; end
        tol = s(1, 1) * numPoles;
        L = rank(s, tol);
    else
        % If numPoles is an integer, use it as the number of poles
        L = round(numPoles);
    end

    % Compute poles using eigenvalue decomposition
    p = eig(inv(s(1:L, 1:L)) * u(:, 1:L)' * A(:, 1:K) * v(:, 1:L));

    % Calculate poles based on time step delta
    delta = t(2) - t(1);
    poles = -log(p) / delta;

    % Filter poles based on the selected method
    n = 1;
    if Method == 'CNRB'
        % CNRB method: keep poles with negative real parts
        for i = 1:length(poles)
            if real(poles(i)) < 0
                pol(n) = poles(i);
                n = n + 1;
            end
        end
    elseif Method == 'CNRU'
        % CNRU method: keep poles with negative real parts and positive imaginary parts
        for i = 1:length(poles)
            if real(poles(i)) < 0 && imag(poles(i)) > 0
                pol(n) = poles(i);
                n = n + 1;
            end
        end
    end
    poles = pol.';  % Return poles as a column vector
end

function [r, x, y] = reconstructMPMint2(poles, t, signal, Method)
    % reconstructMPMint2 reconstructs the signal using the estimated poles
    % and computes the integral of the reconstructed signal.

    N = length(t);  % Length of the time vector
    L = length(poles);  % Number of poles

    % Build the matrix Z based on the poles and time vector
    Z = zeros(N, L);
    for n = 1:N
        for m = 1:L
            Z(n, m) = exp(poles(m) * t(n));
        end
    end

    % Solve the system using least squares
    [rows, cols] = size(signal);
    if rows == 1
        signal = signal';  % Ensure the signal is a column vector
        t = t';
    end
    
   
   r = mldivide(Z, signal);  % CNRB reconstruction method
    
    % Initialize the reconstructed signal and its integral
    x = zeros(N, 1);
    y = zeros(N, 1);
    for i = 1:L
        x = x + ((r(i) * (exp(poles(i) * t))));
        y = y + ((r(i) * (exp(poles(i) * t)) / poles(i)));
    end
    if rows == 1
        x = x';  % Adjust output if original signal was a row vector
        y = y';
    end
    if Method == 'CNRU'
        x = 2*real(x);  % For CNRU method, take the real part
        y = 2*real(y);
    end
end
