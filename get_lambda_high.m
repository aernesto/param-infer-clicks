function L = get_lambda_high(lambda_low, s)
% lambda_low: lambda low
% s: Skellam SNR
% returns lambda high
L=(2 * lambda_low + s^2 + s * sqrt(s^2 + 8 * lambda_low)) / 2;
end
