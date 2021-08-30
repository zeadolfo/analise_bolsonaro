functions {
/* Frank copula log density
*
* @param u Real number on (0,1], not checked but function will return NaN
* @param v Real number on (0,1], not checked but function will return NaN
* @param theta Real number != 0, will throw away otherwise
* @value log density
*/
	real log_frank_copula( real u, real v, real theta) {
		real expr_u = 1.0 - exp( - theta * u );
		real expr_v = 1.0 - exp( - theta * v );
		real expr_theta = 1.0 - exp( - theta );
 // check theta for validity
		if (theta == 0) reject ("theta must be != 0");

		return log( theta * expr_theta * exp(-theta*(u+v)) /
			(expr_theta - (expr_u)*( expr_v) )^2.0 );
	}
}
data {
	int <lower=0> N;
	vector[N] u;
	vector[N] v;
}
parameters {
	real <lower=-1, higher = 1> theta;
}
model {
	vector[N] summands;
	for (n in 1:N) // log -likelihood frank copula density
		summands[n] = log_frank_copula(u[n], v[n], theta );
	target += normal_lpdf( theta | -25.8317 , 10 ); // informative gaussian prior
	target += summands;
}
