##用于实现一些有用的函数

export  add_gaussian_noise, 
        hilbert,
        find_fft_peak

#高斯白噪音注入,信噪比越大越好
function add_gaussian_noise(single::AbstractArray{T, 1}, SNR::Float64) where{T<:Number}
    ps = sum(single.^2) / length(single)
    pn_set = sqrt(ps / (10^(SNR / 10)))

    return single + pn_set * randn(Float64, size(single))
end

function add_gaussian_noise(single::AbstractArray{T, 2}, SNR::Float64; dims=1) where{T<:Number}
    buff = similar(single)
    if dims == 1
        @views Threads.@threads for i in 1:size(single, 2)
            buff[:, i] = add_gaussian_noise(single[:, i], SNR)
        end
    elseif dims == 2
        @views Threads.@threads for i in 1:size(single, 1)
            buff[i, :] = add_gaussian_noise(single[i, :], SNR)
        end
    else
        error("dims should be 1 or 2")
    end
    return @view buff[:,:]
end
#希尔伯特变换实现
function hilbert(x::AbstractArray{T, 1}) where{T<:Real}
    buffer = similar(x, Complex)
    buffer = fft(x, 1)
    len = length(x)
    if len %2 ==0 
        buffer[2: div(len, 2)] .= 2 .*buffer[2: div(len, 2)] 
        buffer[div(len, 2)+1: end] .= 0
    elseif len % 2 == 1 
         buffer[2: div(len+1, 2)] .= 2 .*buffer[2: div(len+1, 2)] 
         buffer[div(len+1, 2)+1:end] .= 0
    end
    buffer = ifft(buffer, 1)
    buffer1 = similar(x, Complex)
    map!((x, y) -> complex(x, imag(y)), buffer1, x, buffer)
    return view(buffer1,:)
end


function hilbert(x::AbstractArray{T, 2}; dims = 1) where{T <: Real}
    buff =  similar(x, Complex)
    if dims == 1
        @views Threads.@threads for i in 1:size(x, 2)
             buff[:, i] = hilbert(x[:, i])
        end
    elseif  dims == 2
        @views Threads.@threads for i in 1:size(x, 1)
             buff[i, :] = hilbert(x[i, :])
        end
    else
        error("dims should be 1 or 2")
    end

    return @view buff[:,:]
end

##进行傅里叶变换并寻找波峰值,找到前5个波峰
##如果调用上面那种阈值需要调整
function find_fft_peak(x::Vector{T}, fs, Threshold::Float64) where{T}
    xf = abs.(fft(x))
    index = Int[]
    L = length(xf)
    P = (xf ./ L)[1:Int(L/2)]
    @. P[2: end] = 2 * P[2: end]
    
    f = view(fs*(0:L-1) / L, 1:Int(L/2))

    ##直接for循环遍历
    Threads.@threads for i in 1:length(P)
        P[i] >= Threshold && push!(index, i)
    end

    return view(f, index)
end
function find_fft_peak(x::Vector,fs, n::Int)
    xf = abs.(fft(x))
    L = length(x)
    index = partialsortperm(view(xf, 1:Int(L/2)), 1:n, rev = true)
    f = (fs*(0:L-1)/L)[1:Int(L/2)]

    return view(f, index)
end