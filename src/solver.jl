##--------------------------------##
##编写求解器函数
##1、SVD-l1
##--------------------------------##
#找到占总比重大于weight前几个特征值,svd分解后s本身已经是从大到小排序完的
export solver

function find_num_max(s, weights)
    @assert typeof(s) <: Vector
    ss = abs.(s)
    weight = 0.0
    total = sum(ss)
    j = 0::Int
    for i in ss
       weight += i/total
       j += 1
       weight >= weights && break
    end
    return j
end

##求解Postion结构体的DOA
function solver(x::Position, fast)
    ex = Expr(:call, x.solver, x, fast)
    eval(ex)
end
##调用前先将data关于x的元素全部清0
function svd_l1!(data::Position,fast = false; weight = 0.95, hyperparameter = 0.8, )

    @assert size(data.real_A, 2) == size(data.real_x, 1)
    @assert size(data.real_A, 2) == size(data.imag_x, 1)
    
    real_A = view(data.real_A, :, :)
    imag_A = view(data.iamg_A, :, :)
    y = data.y

    f = svd(y)
    num_K = find_num_max(f.S, weight)
    IK = zeros(eltype(real_A), num_Ksize(real_A, 1), num_K)
    Threads.@threads for i in 1:num_K
        IK[i, i] = one(eltype(IK))
    end
    Ysv = y* f.V* IK
    real_Ysv = real.(Ysv)
    imag_Ysv = imag.(Ysv)
    n = size(real_A, 2)

    if fast
       ##采用局部细化网格加速计算，但我也不知道有没有问题hh
       ##在原本划分网格的程度上每隔十个点下采样，然后再是原始网格
       seq = Array{Bool, 1}(undef, n)
       fill!(seq, false)
        for i in 1:2
            if i == 1
                real_A1 = view(real_A, :, 1:10:n)
                imag_A1 = view(imag_A, :, 1:10:n)
            elseif i == 2
                real_A1 = view(real_A, :, seq)
                imag_A1 = view(imag_A, :, seq)
            end
            n1 = size(real_A1, 2) 

            model = JuMP.Model(COSMO.Optimizer)
            @variable(model, s_real[1:n1, 1:num_K])
            @variable(model, s_imag[1:n1, 1:num_K])
            @variable(model, 0 <= p)
            @variable(model, z_real[1:length(Ysv)])
            @variable(model, z_imag[1:length(Ysv)])
            @constraint(model, z_real .== (real_Ysv - real_A1*s_real + imag_A1*s_imag)[1:end])
            @constraint(model, z_imag .== (imag_Ysv - imag_A1*s_real - real_A1*s_imag)[1:end])
            @constraint(model, [p; z_real; z_imag] in SecondOrderCone())
            
            @variable(model, r[1:n1])
            for i in 1:n1 #
                @constraint(model, [r[i]; s_real[i, :]; s_imag[i, :]] in SecondOrderCone())
            end
            @variable(model, 0 <= q)
            @constraint(model, q >= sum(r))
            @objective(model, Min, p^2+hyperparameter*q)
            JuMP.optimize!(model)

            if i == 1
                let 
                    Ssv = sum(sqrt.((value.(s_real)).^2 + (value.(s_imag)).^2), dims = 2)
                    power = 10*log10.(Ssv ./ findmax(Ssv)[1]) ##功率大于-10的前几个点的周围
                    index = partialsortperm(power[:, 1], 1:(num_K+1), rev = true) 
                    @views Threads.@threads for i in 1:length(index)
                        max = index[i]+5 > n ? n : index[i] + 5
                        min = index[i]-5 < 1 ? 1 : index[i] - 5
                        @. seq[min:1:max] = true
                    end
                end
            elseif i == 2
                ##data中x的大小肯定大于细化网格结果大小
                ##将计算结果置入data中，重新计算seq
                @. data.real_x[seq, 1:num_K] = value(s_real)
                @. data.imag_x[seq, 1:num_K] = value(s_imag)  
                let
                    Ssv = sum(sqrt.((value.(s_real)).^2 + (value.(s_imag)).^2), dims = 2)
                    power = 10*log10.(Ssv ./ findmax(Ssv)[1]) ##功率大于-10的前几个点的周围
                    index = partialsortperm(power[:, 1], 1:num_K, rev = true)
                    fill!(seq, false)
                    @. seq[index] = true
                end
                        
            else
                error("i should be 1 or 2")
            end

        end
        return seq #返回需要的序列用于绘图
    else
         #那这每次运行模型都是重新，吃完晚饭回来去JuMP文档中查看更新参数的功能
        
         model = JuMP.Model(COSMO.Optimizer)
         @variable(model, s_real[1:n, 1:num_K])
         @variable(model, s_imag[1:n, 1:num_K])
         @variable(model, 0 <= p)
         @variable(model, z_real[1:length(Ysv)])
         @variable(model, z_imag[1:length(Ysv)])
         @constraint(model, z_real .== (real_Ysv - real_A*s_real + imag_A*s_imag)[1:end])
         @constraint(model, z_imag .== (imag_Ysv - imag_A*s_real - real_A*s_imag)[1:end])
 
         @constraint(model, [p; z_real; z_imag] in SecondOrderCone())
 
         @variable(model, r[1:n])
         for i in 1:n
            @constraint(model, [r[i]; s_real[i, :]; s_imag[i, :]] in SecondOrderCone())
         end
 
         @variable(model, 0<=q)
         @constraint(model, q >= sum(r))
         @objective(model, Min, p^2+hyperparameter*q)

         JuMP.optimize!(model)
         
        data.real_x[:, 1:num_K] = value.(s_real)
        data.imag_x[:, 1:num_K] = value.(s_imag) #将结果写入
        #  value.(s_imag)
    end

    return true
end