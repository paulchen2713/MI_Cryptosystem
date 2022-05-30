% 
% Multivariate Cryptography --Matsumoto-Imai(MI) cryptosystem, with input
%   limits that only lower-case English alphabet a~z are allowed, but now
%   the message could be any length, no longer to be restricted to 10.
% 
clear;
clc;
global m n alpha alpha1;
%
% larger m, and n require specific big number calculation, (2^m)^n < 2^52
% 
m = 5;  % GF(2^m)
n = 10;
%
q = 2^m;
one  = 2^m - 1;
zero = 2^m;
%
% generation of qq(x)
%
qq = primitive_polynomial(m);
qqs = size(qq, 1); 
qqs = floor(qqs * rand(1) + 1);
qq = qq(qqs, :);
%
alpha = zeros(2^m, m);
alpha(1, 2) = 1; % alpha^1
for i = 2 : one % i = 2 : 2^m - 1
    if alpha(i-1, m) == 1
        alpha(i, 2:m) = alpha(i-1, 1 : m-1);
        alpha(i, :) = bitxor(alpha(i, :), qq(1:m));
    else
        alpha(i, 2:m) = alpha(i-1, 1 : m-1);
    end
end
%
alphaA = char(alpha + 48);
alphaA = bin2dec(alphaA);
alpha1 = zeros(2^m, 1);
for i = 1 : one % i = 1 : 2^m - 1
    alpha1(alphaA(i) + 1) = i;
end
alpha1(1) = zero;
% 
% generation of g(x)
% 
g = primitive_polynomial(n);
gs = size(g, 1);
gs = floor(gs * rand(1) + 1);
g = g(gs, :);
%
beta = zeros(2^n, n);
beta(1, 2) = 1;
for i = 2 : 2^n - 1
    if beta(i-1, n) == 1
        beta(i, 2:n) = beta(i-1, 1 : n-1);
        beta(i, :) = bitxor(beta(i, :), g(1:n));
    else
        beta(i, 2:n) = beta(i-1, 1 : n-1);
    end
end
% 
% generation of L2, L2I, and cL2
% 
eyee = zero * ones(n, n);
for in = 1 : n
    eyee(in, in) = one;
end
indexL2 = 0; % L2 flag set to 0
while indexL2 == 0
    L2 = floor(zero * rand(n, n) + 1);
    L2I = [L2, eyee];
    L2I = reduced_row_echelon_form_power(L2I);
    L2I = L2I(:, n+1 : 2*n);
    L2_L2I = matrix_multiplication_power(L2, L2I);
    if any(any(L2_L2I - eyee)) == 0
        indexL2 = 1; % L2 flag set to 1
    end
end
cL2 = floor(zero * rand(n, 1) + 1);
%
% generation of L1, L1IM and cL1
%
indexL1 = 0; % L1 flag set to 0
while indexL1 == 0
    L1 = floor(zero * rand(n, n) + 1);
    L1I = [L1, eyee];
    L1I = reduced_row_echelon_form_power(L1I);
    L1I = L1I(:, n+1 : 2*n);
    L1_L1I = matrix_multiplication_power(L1, L1I);
    if any(any(L1_L1I - eyee)) == 0
        indexL1 = 1; % L1 flag set to 1
    end
end
cL1 = floor(zero * rand(n, 1) + 1);
% 
% generation of £c 
% 
ddd = 0;
while ddd ~= 1
    theta = floor((n-1) * rand(1)) + 1;
    tt = q^theta;
    %
    % multiplicative inverse of tt + 1 : t
    %
    [ddd, t] = multiplicative_inverse_p(tt + 1, q^n - 1);
    if ddd ~= 1
        fprintf('Warning for the generation of the multiplicative inverse parameter t !\n\n');
    end
end
% 
% compute LLtt = LL2 ^ tt
% 
LL2 = [L2, cL2];
LLtt = zero * ones(n, n+1);
LLtt(1, :) = LL2(1, :);
for in = 1 : n-1
    int = mod(in * tt, 2^n - 1);
    if int == 0
        bb = beta(2^n - 1, :);
    else
        bb = beta(int, :);
    end
    for ib = 1 : n
        if bb(ib) == 1
            LLtt(ib, :) = matrix_addition_power(LLtt(ib, :), LL2(in + 1, :));
        end
    end
end
%
% compute Ft = LLtt * LL2
%
Ft  = zero * ones(n, n, n);
Ft1 = zero * ones(1, n, n);
Ftc = zero * ones(1, 1, n);
for in1 = 1 : n
    for in2 = 1 : n
        in = mod(in1 + in2 - 2, 2^n - 1);
        if in == 0
            bb = beta(2^n - 1, :);
        else
            bb = beta(in, :);
        end
        for ib = 1 : n
            if bb(ib) == 1
                LLL = matrix_multiplication_power(LLtt(in1, 1:n)', LL2(in2, 1:n));
                Ft(:, :, ib) = matrix_addition_power(Ft(:, :, ib), LLL);
            end
        end
        for ib = 1 : n
            if bb(ib) == 1
                LLL1 = scalar_matrix_multiplication_power(LLtt(in1, n + 1), LL2(in2, 1:n));
                Ft1(1, :, ib) = matrix_addition_power(Ft1(1, :, ib), LLL1);
            end
        end
        for ib = 1 : n
            if bb(ib) == 1
                LLL1 = scalar_matrix_multiplication_power(LL2(in2, n + 1), LLtt(in1, 1:n));
                Ft1(1, :, ib) = matrix_addition_power(Ft1(1, :, ib), LLL1);
            end
        end
        mm = multiplication_power(LLtt(in1, n + 1), LL2(in2, n + 1));
        for ib = 1 : n
            if bb(ib) == 1
                Ftc(1, 1, ib) = addition_power(Ftc(1, 1, ib), mm);
            end
        end
    end
end
%
% compute FB = L1(Ft)
%
FB  = zero * ones(n, n, n);
FB1 = zero * ones(1, n, n);
FBc = zero * ones(1, 1, n);
for in = 1 : n
    for inn = 1 : n
        LF = scalar_matrix_multiplication_power(L1(in, inn), Ft(:, :, inn));
        FB(:, :, in) = matrix_addition_power(FB(:, :, in), LF);
    end
    for inn = 1 : n
        LF1 = scalar_matrix_multiplication_power(L1(in, inn), Ft1(1, :, inn)); %
        FB1(1, :, in) = matrix_addition_power(FB1(1, :, in), LF1);
    end
    LFc(1 : n, 1) = Ftc(1, 1, :);
    LFc = matrix_multiplication_power(L1(in, :), LFc);
    FBc(1, 1, in) = matrix_addition_power(FBc(1, 1, in), LFc);
    FBc(1, 1, in) = addition_power(FBc(1, 1, in), cL1(in));
end
%
% public key:
%     m, n, FB, FB1, FBc
% private key:
%     theta, L1, cL1, L2, cL2
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
% Encryption
%
% in this case, the message input could be any length, but still only 
% lower-case English alphabet allowed. no upper-case, no Greek alphabet, 
% no digit 0~9, no spaces. surprisingly, '~', '{}', '|'symbol works fine.
%
% testing:
%   '~!@#$%^&*()_+[]{}\|;''",./<>? 0123456789 QWERTYUIOPASDFGHJKLZXCVBNM', 
%   'qwertyuiopasdfghjklzxcvbnm'
%   'not~exactly~ten~letters~is~allowed~now~but~still~no~uppercase~no~spaces~no~digits'
% 
message = 'not~exactly~ten~letters~is~allowed~now~but~still~no~uppercase~no~spaces~no~digits';
message_len = length(message);
mN = ceil(message_len / n);
mR = mod(message_len, n);
num_a = 0;
if mR ~= 0
    for i = mR + 1 : n
        message = strcat(message, 'a');
        num_a = num_a + 1;
    end
end
xx = double(message) - 96;
%
% compute the cipher = FB(x)
%
cipher = zero * ones(n, mN);
for im = 1 : mN
    x = (xx((im - 1)*n + 1 : im * n))';
    for in = 1 : n
        cipher(in, im) = addition_power(cipher(in, im), ...
            matrix_multiplication_power(matrix_multiplication_power(x', FB(:, :, in)), x));
        cipher(in, im) = addition_power(cipher(in, im), ...
            matrix_multiplication_power(FB1(1, :, in), x));
        cipher(in, im) = addition_power(cipher(in, im), FBc(in));
    end
end
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
% Decryption
%
yy = cipher;
%
% express t = q^t1 + q^t2 + .. + q^tk + tr, q = 2^m
%
T = zeros(1, 0);
index = 0;
while t >= q
    index = index + 1;
    ii = 0;
    while q^(ii + 1) <= t
        ii = ii + 1;
    end
    T(index) = q^ii;
    t = t - q^ii;
end
if t ~= 0
    T(index + 1) = t;
end
%
% compute L1y = L1^(-1)(y), inverse_L1(y)
%
mN = size(yy, 2);
message_R = char();
for im = 1 : mN
    y = yy(:, im);
    y = matrix_addition_power(y, cL1);
    zB = matrix_multiplication_power(L1I, y);
    %
    %
    Ty = size(T, 2); % set Ty = k + 1
    FtI = zero * ones(n, Ty);
    for it = 1 : Ty-1 %
        at = zero * ones(n, 1);
        at(1) = zB(1);
        for in = 1 : n-1
            int = mod(in * T(it), 2^n - 1);
            if int == 0
                bb = beta(2^n - 1, :);
            else
                bb = beta(int, :);
            end
            for ib = 1 : n
                if bb(ib) == 1 %
                    at(ib) = addition_power(at(ib), zB(in + 1));
                end
            end
        end
        FtI(:, it) = at;
    end
    %
    if mod(T(Ty), zero) == 0
        at = zero * ones(n, 1);
        at(1) = zB(1);
        for in = 1 : n-1
            int = mod(in * T(Ty), 2^n - 1);
            if int == 0
                bb = beta(2^n - 1, :);
            else
                bb = beta(int, :);
            end
            for ib = 1 : n
                if bb(ib) == 1 %
                    at(ib) = addition_power(at(ib), zB(in + 1));
                end
            end
        end
        FtI(:, Ty) = at;
    else
        at = L1y;
        for it = 1 : T(Ty)-1
            temp = at;
            at = zero * ones(n, 1);
            for in1 = 0 : n-1
                for in2 = 0 : n-1
                    in = mod(in1 + in2, 2^n - 1);
                    if in == 0
                        bb = beta(2^n - 1, :);
                    else
                        bb = beta(in, :);
                    end
                    mm = multiplication_power(temp(in1 + 1), zB(in2 + 1));
                    for ib = 1 : n
                        if bb(ib) == 1
                            at(ib) = addition_power(at(ib), mm);
                        end
                    end
                end
            end
        end
        FtI(:, Ty) = at;
    end
    %
    % compute z = Ft^(-1)(L1y)
    %
    z = FtI(:, 1);
    for it = 1 : Ty-1
        temp = z;
        z = zero * ones(n, 1);
        for in1 = 0 : n-1
            for in2 = 0 : n-1
                in = mod(in1 + in2, 2^n - 1);
                if in == 0
                    bb = beta(2^n - 1, :);
                else
                    bb = beta(in, :);
                end
                mm = multiplication_power(temp(in1 + 1), FtI(in2 + 1, it + 1));
                for ib = 1 : n
                    if bb(ib) == 1
                        z(ib) = addition_power(z(ib), mm);
                    end
                end
            end
        end
    end
    %
    % compute x_R = L2^(-1)(z), inverse_L2(z)
    %
    z = matrix_addition_power(z, cL2);
    x_R = matrix_multiplication_power(L2I, z);
    message_R = strcat(message_R, (char(x_R + 96))');
end
%
% R_len = length(message_R);
% message_R = message_R(1 : R_len-num_a);
%
fprintf('the original message with concatenated ''a'' is: %s\n', message);
fprintf('the recovery message with concatenated ''a'' is: %s\n\n', message_R);
fprintf('the original message is: %s\n', message(1 : message_len));
fprintf('the recovery message is: %s\n\n', message_R(1 : message_len));
%
%


