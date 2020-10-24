%% test
fun = @(x) (1/100 * exp(-8*sqrt(3) * x));
integral(fun, 0, 100)