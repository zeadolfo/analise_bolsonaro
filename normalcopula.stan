functions {
/* Frank copula log density
*
* @param u Real number on (0,1], not checked but function will return NaN
* @param v Real number on (0,1], not checked but function will return NaN
* @value log density
*/
	real log_normal_copula( real u, real v, real rho) {
		return log(1/((1 - rho^2)^(-1/2)) * exp((u^2 + v^2 + 2 * rho * u * v)/(2 * (1 - rho^2))) * exp((u^2 + v^2)/2)^2);
	}
	real extended_beta(real rho, real a, real b){
		return log((1 + rho)^(a - 1) * (1 - rho)^(b - 1));
	}
}
data {
	int <lower=0> N;
	vector[N] u;
	vector[N] v;
}
parameters {
	real <lower=-1, upper = 1> rho;
}
model {
	vector[N] summands;
	for (n in 1:N) // log -likelihood frank copula density
		summands[n] = log_normal_copula(u[n], v[n], rho);
	target += extended_beta(rho, 1 , 1 ); // non-informative gaussian prior
	target += summands;
}
