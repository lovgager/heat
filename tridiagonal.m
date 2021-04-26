function x = tridiagonal(a, b, d)
%TRIDIAGONAL a(double) - upper & lower diags, b(double) - main diag,
%d(array of double) - right part of the system
    dim = size(d);
    dim = dim(2);
    alpha = zeros(1, dim); alpha(1) = -a/b;
    beta = zeros(1, dim); beta(1) = d(1)/b;
    gamma = zeros(1, dim); gamma(1) = b;

    for i = 2:dim
        gamma(i) = b + a * alpha(i - 1);
        beta(i) = (d(i) - a * beta(i - 1)) / gamma(i);
        alpha(i) = -a / gamma(i);
    end

    x = zeros(1, dim);
    x(dim) = beta(dim);
    for i = dim - 1:-1:1
        x(i) = alpha(i) * x(i + 1) + beta(i);
    end
end

