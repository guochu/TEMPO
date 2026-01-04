
function pt_ti_mpotensor(corr::CorrelationMatrix, op1::AbstractMatrix, op2::AbstractMatrix, alg::ExponentialExpansionAlgorithm)
	m1 = GenericDecayTerm(op1, op2, corr.ηₖⱼ[2:end])
	m2 = GenericDecayTerm(op2, op1, corr.ηⱼₖ[2:end])

	m1s = exponential_expansion(m1, alg=alg)
	m2s = exponential_expansion(m2, alg=alg)

	# println("here---", corr.ηₖⱼ[1], " ", corr.ηⱼₖ[1])
	eta = corr.ηₖⱼ[1] + corr.ηⱼₖ[1]
	# h1 = corr.ηₖⱼ[1] * op1 * op2 + corr.ηⱼₖ[1] * op2 * op1
	# h1 = corr.ηₖⱼ[1] * op2 * op1 + corr.ηⱼₖ[1] * op1 * op2

	# h1 = (corr.ηₖⱼ[1] + corr.ηⱼₖ[1])  * op1 * op2 
	h1 = (eta/2) .* (op1 * op2 + op2 * op1)
	return SchurMPOTensor(h1, vcat(m1s, m2s))
	# return SchurMPOTensor(h1, [])
end
