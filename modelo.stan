functions { 
	real log_normal_copula( real u, real v, real rho) {
		return log(1/(sqrt(1 - rho^2)) * exp(-(u^2 + v^2 - 2 * rho * u * v)/(2 * (1 - rho^2))) * exp((u^2 + v^2)/2));
	}
	real extended_beta(real rho, real a, real b) {  
		return log((1+rho)^(a - 1) * (1 - rho)^(b - 1));
	}
}
data {
	int <lower=0> N;
	vector[N] u;
	vector[N] v;
}
parameters {
	real <lower= - 1,upper=1> rho;
}
model {
	vector[N] summands;
	for (n in 1:N)	{
		summands[n] = log_normal_copula(u[n], v[n], rho);
	}
	target += extended_beta(rho, 1, 1); 
	target += summands;
}
