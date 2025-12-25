


function spin_half_matrices()
	s_SP = Array{Float64, 2}([0 0; 1 0])
	s_SM = Array{Float64, 2}([0 1; 0 0])
	s_Z = Array{Float64, 2}([-1 0; 0 1])
	s_x = s_SP+s_SM
	s_y = -im*(s_SP-s_SM)
	n = Array{Float64, 2}([0 0; 0 1])
	return Dict("x"=>s_x, "y"=>s_y, "z"=>s_Z, "+"=>s_SP, "-"=>s_SM, "n"=>n)
end

pauli_x() = Matrix{Float64}([0 1; 1 0])
pauli_y() = Matrix{ComplexF64}([0 im; -im 0])
pauli_z() = Matrix{Float64}([-1 0; 0 1])
spin_up() = Matrix{Float64}([0 0; 0 1])
spin_down() = Matrix{Float64}([1 0; 0 0])


function rabi_ham(Ω; d)
	x = pauli_x()
	z = pauli_z()
	# sp = Array{Float64, 2}([0 0; 1 0])
	# ns = [0 0; 0 1.]
	Is = one(x)
	a = bosonaoperator(d=d)
	adag = a'
	n = bosondensityoperator(d=d)
	Ib = one(n)

	Himp = Ω * kron(x, Ib)
	Hbath = kron(Is, n)
	Hhyb = kron(z, adag+a)

	H = Himp + Hhyb + Hbath
	return H, n
end

function rabi_ham_2(Ω; d)
	x = pauli_x()
	z = pauli_z()
	# sp = Array{Float64, 2}([0 0; 1 0])
	# ns = [0 0; 0 1.]
	Is = one(x)
	a = bosonaoperator(d=d)
	adag = a'
	n = bosondensityoperator(d=d)
	Ib = one(n)

	Himp = Ω * kron(z, Ib)
	Hbath = kron(Is, n)
	Hhyb = kron(x, adag+a)

	H = Himp + Hhyb + Hbath
	return H, n
end

function _rand_dm(d)
	x = randn(ComplexF64, d, d)
	x = x' * x
	return x / tr(x)
end
