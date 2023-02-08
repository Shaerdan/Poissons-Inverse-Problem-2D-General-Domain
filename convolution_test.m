f = [1,2,3,4,5];
kernel = [0.25, 0.5, 0.25; 0.1, 0.1, 0.1];

% Reshape the source vector into a 2D matrix
f_matrix = reshape(f, [length(f), 1]);

% Perform the 2D convolution on the reshaped matrix
f_matrix_conv = conv2(f_matrix, kernel, 'same');

% Reshape the result back into a 1D vector
f_conv = reshape(f_matrix_conv, [length(f), 1]);

%This will perform a 1D convolution in the first dimension (rows) of the 2D matrix. If you want to perform the convolution in the second dimension (columns), you can transpose the source matrix before the convolution and transpose the result back after the convolution.
