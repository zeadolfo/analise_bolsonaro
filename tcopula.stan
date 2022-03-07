functions { 
	real log_t_copula( real u, real v, real rho, real d) {
		return log(((1 + ((u^(2) + v^(2) - 2 * rho * u^(1) * v^(1))/(d * (1 - rho^2))))^(-(d + 2)/2))/(d*(1 - rho^2)/((sqrt(1 - rho^2) * ((1 + ((u^(2))/d))^(-(d + 1)/2)) * (1 + ((v^(2))/d))^(-(d + 1)/2)))))	; 
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
	real <lower=-1,upper=1> rho;
}
model {
	vector[N] summands;
	for (n in 1:N)	{
		summands[n] = log_t_copula(u[n], v[n], rho, 4);
	}
	target += extended_beta(rho, 1, 1); 
	target += summands;
}
