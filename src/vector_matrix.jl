##-------------------------------------------##
##先定义形状喝流行阵列参数,然后在调用函数生成对应A
##-------------------------------------------##
##name = "线路输入 (LIS_MIIS_8X8_192K)"
export linear1,
        linear2,
        circle,
        A_matrix

abstract type ArrayShape end

struct linear1{T1<:Float64, T2<:Int} <: ArrayShape #一维直线阵列
    d::T1  #阵列间距
    n::T2  #阵列个数
end
struct linear2{T1<:Float64, T2<:Int} <: ArrayShape   #二维阵列
    d1::T1
    d2::T1 #阵列间距
    n1::T2 #长度阵列个数
    n2::T2 #宽度阵列个数
end
struct circle{T1<:Float64, T2<:Int} <:ArrayShape #圆形阵列
    n::T2  #阵列个数
    r::T1 #阵列半径
end
struct A_matrix{N,T<:Float64,F<:ArrayShape,T1<:NTuple{N, Float64},T2<:NTuple{N, Float64},I<:Int }
    angle::T1  #角度范围,单位为°
    length::T2  #声源波长,初始化是用Vector声明，不然无法修改数值
    step1::T    #角度步进值1,单位为弧度制
    shape::F    #排列形状
    dims::I     #空间角度划分维度
end

function (A::A_matrix)()
    if typeof(A.shape) == linear1
        result = linear1_methon(A)
    elseif typeof(A.shape) == linear2
        result = linear2_methon(A)
    elseif typeof(A.shape) == circle
        result = circle_methon(A)
    else
        error("shape shouldn't be $(typeof(A.shape))")
    end
    return result
end

#一维直线阵列只能用于测量一维空间角度
function linear1_methon(A::A_matrix)
    angle, bolength, step1, shape = (A.angle .*pi ./ 180), A.length[1], A.step1, A.shape
    d, n = shape.d, shape.n
    @assert A.dims == 1

    d = 0: d: (n-1)*d
    theta = range(angle[1], stop=angle[2], step=step1)
    Am = Array{Complex, 2}(undef, n, length(theta))
    @. Am = exp(-im*2*pi/bolength*d*theta')
    return @view Am[:,:] #返回结果
end

#适用于空间网格二维划分，并且对应空间二维ULA阵列形式
function linear2_methon(A::A_matrix)
    angle, bolength, step1, shape = (A.angle .*pi ./ 180), A.length[1], A.step1, A.shape
    d1, d2, n1, n2 = shape.d1, shape.d2, shape.n1, shape.n2

    @assert length(angle) == 4
    @assert length(step1) == 2
    
    theta = range(angle[1], step=step1[1], stop=angle[2])
    beta = range(angle[3], step=step1[2], stop=angle[4])
    l1, l2 = length(theta), length(beta)
    array_dim1 = 0: d1: (n1-1)*d1
    array_dim2 = 0: d2: (n2-1)*d2
    Am = Array{Complex,2}(undef, n1*n2, l1*l2)
    #生成d的序列,具体X，Y，Z轴的赋值与麦克风安装位置有关
    d_seq_buff = Array{Float64,2}(undef, n1*n2, 3)
    @views Threads.@threads for i in 1:length(array_dim1)
        @. d_seq_buff[(i-1)*n1+1 : (i-1)*n1+n2, 1] = array_dim1[i]
        d_seq_buff[(i-1)*n1+1 : (i-1)*n1+n2, 3] = array_dim2[1:end]
    end
    @. d_seq_buff[:, 2] = 0 
    #生成theta序列
    Trigonometric_buff1 = Array{Float64, 2}(undef, 2, l1) #存储三角函数值
    @. Trigonometric_buff1[1, :] = sin(theta)
    @. Trigonometric_buff1[2, :] = cos(theta)
    Trigonometric_buff2 = Array{Float64, 2}(undef, 2, l2) #存储三角函数值
    @. Trigonometric_buff2[1, :] = sin(beta)
    @. Trigonometric_buff2[2, :] = cos(beta)
    #序列相乘得到Am
    @views Threads.@threads for i in 1: l1
        @. Am[:, (i-1)*l1+1:(i-1)*l1+l2] = (d_seq_buff[:, 1] * Trigonometric_buff2[2, :]' * Trigonometric_buff1[1, i])
                        +(d_seq_buff[:, 3] * Trigonometric_buff1[2, i])
    end
    
    @. Am = exp(-im * 2 * pi / bolength * Am)

    return @view Am[:, :]
end

#圆形阵列，对应空间二维划分
function circle_methon(A::A_matrix)
    angle, bolength, step1, shape = (A.angle .*pi ./ 180), A.length[1], A.step1, A.shape
    n, r = shape.n, shape.r
    @assert length(angle) == 4
    @assert length(step1) == 2
    
    theta = range(angle[1], step=step1[1], stop=angle[2])
    beta = range(angle[3], step=step1[2], stop=angle[4])
    l1, l2 = length(theta), length(beta)
    Am = Array{Complex, 2}(undef, n, l1*l2)
    #获得阵列位置序列
    d_seq_buff = Array{Float64, 2}(undef, n, 3)
    P_buff = (0:n-1) * 2 * pi / n
    @. d_seq_buff[1:n, 1] = sin(P_buff) * r 
    @. d_seq_buff[1:n, 2] = 0
    @. d_seq_buff[1:n, 3] = cos(P_buff) * r
    #获得theta序列
    Trigonometric_buff1 = Array{Float64, 2}(undef, 2, l1) #存储三角函数值
    @. Trigonometric_buff1[1, :] = sin(theta)
    @. Trigonometric_buff1[2, :] = cos(theta)
    Trigonometric_buff2 = Array{Float64, 2}(undef, 2, l2) #存储三角函数值
    @. Trigonometric_buff2[1, :] = sin(beta)
    @. Trigonometric_buff2[2, :] = cos(beta)
    #序列相乘得到Am
    @views Threads.@threads for i in 1: l1
        @. Am[:, (i-1)*l1+1:(i-1)*l1+l2] = (d_seq_buff[:, 1] * Trigonometric_buff2[2, :]' * Trigonometric_buff1[1, i])
                        +(d_seq_buff[:, 3] * Trigonometric_buff1[2, i])
    end
    @. Am = exp(-im * 2 * pi / bolength * Am)

    return @view Am[:, :]
end

