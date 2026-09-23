

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

# ---- ⟨B| m |A⟩ 三链环境原语（PT：rank-4 site 张量）----
# 与 FiniteMPSAlgorithms mult 引擎的三链环境约定一致（hstorage 轴序 (B, m, A)），
# 供 mult 的 finalize 以及 TDVPIF 的 PT 流使用。

function updateleft(cleft::DenseMPSTensor, B::DenseMPOTensor, m::DenseMPOTensor, A::DenseMPOTensor)
    @tensor tmp[9,8,5] := ((cleft[1,2,3] * A[3,4,5,6]) * m[2,7,8,4]) * conj(B[1,7,9,6])
    return tmp
end

function updateright(cright::DenseMPSTensor, B::DenseMPOTensor, m::DenseMPOTensor, A::DenseMPOTensor)
    @tensor tmp[1,7,9] := ((conj(B[1,2,3,4]) * cright[3,5,6]) * m[7,2,5,8] ) * A[9,8,6,4]
    return tmp
end

function reduceH_single_site(A::DenseMPOTensor, m::DenseMPOTensor, cleft::DenseMPSTensor, cright::DenseMPSTensor)
	@tensor tmp[1,7,9,6] := ((cleft[1,2,3] * A[3,4,5,6]) * m[2,7,8,4]) * cright[9,8,5]
    return tmp
end
