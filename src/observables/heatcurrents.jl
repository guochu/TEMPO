

function heatcurrent(lattice::RealADTLattice1Order)
	
end


function build_current_op(lattice::RealADTLattice1Order, corr::RealCorrelationFunction, hyb::AdditiveHyb, k::Int)
	row = index(lattice, k, branch=:+)
	cols = Int[]
	T = promote_type(scalartype(lattice), scalartype(corr))
	coefs = T[]
	for j in k-1:-1:1
		for  b2 in branches(lattice)
			col = index(lattice, j, branch=b2)
			push!(cols, col1)
			push!(coefs, index(corr, k, j, b1=:+, b2=b2))
		end
	end
end


function build_current_op(lattice::RealPTLattice1Order, corr::RealCorrelationFunction, hyb::NonAdditiveHyb, k::Int)
	
end