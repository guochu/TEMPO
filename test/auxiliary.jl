println("------------------------------------")
println("|            Auxiliary             |")
println("------------------------------------")


@testset "Exponential expansion    " begin
	L = 100
	atol = 1.0e-5
	for alpha in (-2, -2.5, -3)
		xdata = [convert(Float64, i) for i in 1:L]
		ydata = [1.3 * x^alpha for x in xdata]
		xs1, lambdas1 = exponential_expansion(ydata, PronyExpansion(n=20,tol=atol))
		@test expansion_error(ydata, xs1, lambdas1) < atol
		xs2, lambdas2 = exponential_expansion(ydata, DeterminedPronyExpansion(n=20,tol=atol))
		@test expansion_error(ydata, xs2, lambdas2) < atol
		# xs3, lambdas3 = exponential_expansion(ydata, PronyExpansion2(n=20,atol=atol))
		# @test expansion_error(ydata, xs3, lambdas3) < atol
		# xs2, lambdas2 = exponential_expansion(ydata, LsqExpansion(atol=atol))
		# @test expansion_error(ydata, xs2, lambdas2) < atol
	end
	L = 500
	xdata = [convert(Float64, i) for i in 1:L]
	ydata = [1.3 * 0.7^x + 0.7 * 0.5^x - 1.1 * 0.8^x + 0.1*0.95^x for x in xdata]
	xs1, lambdas1 = exponential_expansion(ydata, PronyExpansion(n=20,tol=atol))
	@test expansion_error(ydata, xs1, lambdas1) < atol
	for stepsize in 2:6
		xs2, lambdas2 = exponential_expansion(ydata[1:stepsize:L], PronyExpansion(n=20, stepsize=stepsize, tol=atol))
		@test expansion_error(ydata, xs2, lambdas2) < atol
	end
	L = 500
	xdata = [convert(Float64, i) for i in 1:L]
	ydata = [(1.3+0.2im) * (0.7+0.3im)^x + (0.7+1.1im) * (0.5+0.1im)^x - (1.1-0.3im) * (0.8-0.5im)^x + (0.1-0.2)*(0.95+0.2im)^x for x in xdata]
	for alg in (PronyExpansion(n=20,tol=atol),)
		xs, lambdas = exponential_expansion(ydata, alg)
		@test expansion_error(ydata, xs, lambdas) < atol
	end
end