

"""
	updateright(hold::AbstractArray{T, 2}, mpsAj::MPSTensor, mpsBj::MPSTensor) where {T<:Number}
update storage from right to left for overlap of mps
"""
function updateright(hold::AbstractMatrix, mpsAj::DenseMPSTensor, mpsBj::DenseMPSTensor)
	@tensor m2[-1 -2;-3] := conj(mpsAj[-1, -2, 1]) * hold[1, -3]
	@tensor hnew[-1;-2] := m2[-1,1,2] * mpsBj[-2,1,2]
	return hnew
end

"""
	updateleft(hold::AbstractArray{T, 2}, mpsAj::MPSTensor, mpsBj::MPSTensor) where {T<:Number}
update storage from left to right for overlap of mps
"""
function updateleft(hold::AbstractMatrix, mpsAj::DenseMPSTensor, mpsBj::DenseMPSTensor)
	@tensor m2[-1 -2; -3] := conj(mpsAj[1, -2, -3]) * hold[1, -1]
	@tensor hnew[-1; -2] := m2[1,2,-1] * mpsBj[1,2,-2]
	return hnew
end


function n_fuse(m::AbstractArray{<:Number, N}, i::Int) where {N}
	@assert (1 <= i < N)
	@assert size(m, i) == size(m, i+1)
	# @assert (i < M) || (M <= i < N)
	local tmp
	s_front = ntuple(x->size(m, x), i-1)
	s_tail = ntuple(x->size(m, i+1+x), N-i-1)
	s1 = prod(s_front)
	s2 = prod(s_tail)

	d = size(m, i)
	m4 = reshape(m, s1, d, d, s2)
	m3 = zeros(scalartype(m), s1, d, s2)
	# m3[:,1,:] = m4[:, 1, 1, :]
	# m3[:,2,:] = m4[:, 1, 2, :] + m4[:, 2, 1, :] + m4[:, 2, 2, :]

	for j in 1:d
		m3[:,j,:] = m4[:, j, j, :]
	end


	new_size = (s_front..., d, s_tail...)
	return reshape(m3, new_size)
end

# ---- ⟨z| A ⊙ B⟩ 三链环境原语（ADT：物理指标共享的逐点乘积）----
# 与 FiniteMPSAlgorithms hadamard 引擎的环境约定一致（hstorage 轴序 (z, A, B)），
# 供 mult 的 finalize 以及 TDVPIF 的 ADT 流使用。

function updatemultleft(left::DenseMPSTensor, zj::DenseMPSTensor, xj::DenseMPSTensor, yj::DenseMPSTensor)
	@tensor tmp1[1,5,4;2] := left[1,2,3] * yj[3,4,5]
	@tensor tmp2[1,3,5;6,2] := tmp1[1,2,3,4] * xj[4,5,6]
	tmp3 = n_fuse(tmp2, 2)
	@tensor tmp2[5,3;4] := tmp3[1,2,3,4] * conj(zj[1,2,5])
	return tmp2
end

function updatemultright(right::DenseMPSTensor, zj::DenseMPSTensor, xj::DenseMPSTensor, yj::DenseMPSTensor)
	@tensor tmp1[4; 1 2 5] := yj[1,2,3] * right[3,4,5]
	@tensor tmp2[4 1 2 5; 6] := xj[1,2,3] * tmp1[3,4,5,6]
	tmp3 = n_fuse(tmp2, 3)
	@tensor tmp2[4,5;1] := conj(zj[1,2,3]) * tmp3[4,5,2,3]
	return tmp2
end


function get_left_xy(left::DenseMPSTensor, xj::DenseMPSTensor, yj::DenseMPSTensor)
	@tensor tmp1[1,5,4;2] := left[1,2,3] * yj[3,4,5]
	@tensor tmp2[1,3,5;6,2] := tmp1[1,2,3,4] * xj[4,5,6]
	tmp3 = n_fuse(tmp2, 2)
	return tmp3
end
function get_xy_right(right::DenseMPSTensor, xj::DenseMPSTensor, yj::DenseMPSTensor)
	@tensor tmp1[4; 1 2 5] := yj[1,2,3] * right[3,4,5]
	@tensor tmp2[4 1 2 5; 6] := xj[1,2,3] * tmp1[3,4,5,6]
	tmp3 = n_fuse(tmp2, 3)
	return tmp3
end
