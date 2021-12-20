##此文件用来进行定位解算的
##先写下测试文件，保证A，y和求解器能正常工作
export Position

struct Position{N, N1<:Int, F<:Symbol, T1<:Array{Float64, N}, T2<:Array{Float64, 2}, T3<:Array{Complex, 2}}
    y::T3
    real_x::T1    #求得的解
    imag_x::T1
    real_A::T2    #过完备原子基
    iamg_A::T2
    dims::N1   #定位维度
    solver::F    #求解方法，采用元编程，便于后期扩展
end
#角度方位，空间网格维度和步进角度,中心频率和波速度
function Position(y, A, dims::Int, solver::Symbol) where{N}
    @assert size(y, 1) == size(A, 1)
    real_x = Array{Float64, 2}(undef, size(A, 2), size(y, 2))
    imag_x = Array{Float64, 2}(undef, size(A, 2), size(y, 2))
    fill!(real_x, zero(eltype(real_x)))
    fill!(imag_x, zero(eltype(imag_x)))
    return Position(y, real_x, imag_x, real.(A), imag.(A), dims, solver)
end



