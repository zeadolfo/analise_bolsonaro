functions {
	real log_frank_copula( real u, real v, real theta) {
		real expr_u = 1.0 - exp( - theta * u );
		real expr_v = 1.0 - exp( - theta * v );
		real expr_theta = 1.0 - exp( - theta );
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
	real <lower=-100, upper = 100> theta;
}

model {
	vector[N] summands;
	for (n in 1:N) 
		summands[n] = log_frank_copula(u[n], v[n], theta );
	target += normal_lpdf( theta | -100 , 10		); 
	target += summands;
}
