

function updateright(hold::AbstractMatrix, hAj::DenseMPOTensor, hBj::DenseMPOTensor) 
	@tensor m2[-1 -2;-3 -4] := conj(hAj[-1, -2, 1, -4]) * hold[1, -3]
	@tensor hnew[-1;-2] := m2[-1,1,2,3] * hBj[-2,1,2,3]
	return hnew
end

function updateleft(hold::AbstractMatrix, hAj::DenseMPOTensor, hBj::DenseMPOTensor) 
	@tensor m2[-1 -2 ; -3 -4] := conj(hAj[1, -2, -3, -4]) * hold[1, -1]
	@tensor hnew[-1; -2] := m2[1,2,-1,3] * hBj[1,2,-2,3]
	return hnew
end