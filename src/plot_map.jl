##绘制定位图像,其实就是左边是一个绘制定位图像，右边是绘制时域图像
##定位图像每帧都是新的，时域图像则是在之前的基础上进行修改
##for循环中直接plot将不会展示图片，Plots中图片只有在return中才会被展示
##想要在plot中直接展示需要添加show = true
export  plot_map,
        map_plot,
        pc, sp
##两幅图，一个是定位图，一个是实时接受的图,使用flag判断是否画图
struct map_plot{T1<:Vector{Plots.Plot}, T2<:Vector{Bool}, T3<:Float64}
    pc::T1
    flag::T2
    t::T3
end

##定义画板的一些参数并配置对应的修改函数
Base.@kwdef mutable struct plot_parameters
    size::NTuple{Int, 2} = (800, 600)
end
parameter = plot_parameters()
#自定义画板大小
function set_size(L, W)
    @assert typeof(L) == Int
    @assert typeof(W) == Int
    parameter.size = (L,W)
end


##定义整个参数类型，然后调用default对画板进行初始化
default(size=parameter.size)

function map_plot(;N=2)
    map_plot(Vector{Plots.Plot}(undef, N), Vector{Bool}(undef, N), 0.0)
end
#是否只更新定位结果
function map_plot_init(n::Int;show_time = true)
    @assert n>= 2

    pct = map_plot(;N=n)
    fill!(pct.flag, true)
    show_time == true ? nothing : pct.flag[2] = false
    for p in pct.pc
        p = Plots.plot() #不加直接layout回报错
    end
    return pct
end
#直接生成两个figure用于绘图，生成后就能调用heatmap!和plot!函数不报错
const pc = map_plot_init(2)
const sp = Plots.plot(pc.pc..., layout=(1,2))
Plots.plot!(sp[1], annotations = (0.5, 0.5, Plots.text("Position Map", :center)))
Plots.plot!(sp[2], annotations = (0.5, 0.5, Plots.text("Collect data", :center)))
Plots.plot!(sp, size=(800, 600))
##是否展示收集到的时域信号
##Vector => Matrix, 返回很横坐标和纵坐标
function vector2matrix(x::Vector{Float64}, a::A_matrix)
    theta = range(a.angle[1], step=a.step1[1], a.angle[2])
    beta = range(a.angle[3], step=a.step1[1], a.angle[4])
    l1, l2 = length(theta), length(beta)

    power = Array{Float64, 2}(undef, l2, l1)
    @. power = reshape(x, l2, l1)
    return view(power,:,:) , (theta, beta)
end

function plot_map(x::Position, vm::A_matrix; time_img = false)
    if time_img
        plot_time_signal(x.y)
    end
    #绘制定位图像用heatmap,x的power表示图像，因为这个有可能为0，就不采用对数形式
    #先将复数形式转换为实数形式
    power_buff = sum(sqrt.(x.imag_x.^2 + x.imag_y.^2), dims=2)
    if vm.dims == 1
        angle = range(vm.angle[1], step=vm.step1[1], step=vm.angle[2])
        Plots.plot(sp[1], angle, power_buff)
    elseif vm.dims == 2
        power, angle = vector2matrix(power_buff, vm)
        Plots.heatmap!(sp[1], angle[1], angle[2], power)
    else
        println("只进行一维或者二维的绘图")
        return false
    end
    return true
end

##绘制复数图像
function plot_time_signal(x::Array{Complex, N}) where{N}
    plot_time_signal(real.(x))
end

function plot_time_signal(x::Array{T, N}) where{T<:Real, N}
    @assert eltype(x) <: Real
    plot(sp[2], 1:length(x), x, lab="real data")
end
