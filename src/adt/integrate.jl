
function integrate(x::ADT)
	L = length(x)
	sca = scaling(x)
	v = dropdims(sum(x[L], dims=2), dims=(2,3)) * sca
	for i in L-1:-1:1
		tmp = dropdims(sum(x[i], dims=2), dims=2) * sca
		v = tmp * v
	end
	return only(v)
end


function integrate(x::ADT, y::ADT)
	(length(x) == length(y)) || throw(DimensionMismatch("adt size mismatch"))
	sca = scaling(x) * scaling(y)
	L = length(x)
	@tensor v[1,4] := sca * x[L][1,2,3] * y[L][4,2,3]
	for i in L-1:-1:1
		@tensor tmp[1,4] := sca * x[i][1,2,3] * y[i][4,2,5] * v[3,5]
		v = tmp
	end
	return tr(v)
end