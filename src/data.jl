##------------------------------##
##实现从下位机读取麦克风阵列数据
##1、调用MATLAB函数读取
##2、c++封装的SO读取
##3、自己生成随机数
##------------------------------##
##name = "线路输入 (LIS_MIIS_8X8_192K)"
export  MicrophoneParameter,
        init_matlab,
        getdata

struct MicrophoneParameter{T1<:String, T2<:Int} 
    name::T1 #麦克风阵列名称
    fs::T2  #采样频率
    sample::T2  #一帧采样数
    numchannels::T2    #通道数
end

function init_matlab(parameter::MicrophoneParameter)
    s = MATLAB.MSession(0)
    put_variable(s, :name, parameter.name)
    eval_string(s, "H = audioDeviceReader")
    eval_string(s, "H.Device = name ")
    eval_string(s, "H.NumChannels = $(parameter.numchannels)") #通道数要跟实际设备中搜集到的数据一致
    eval_string(s, "H.SampleRate = $(parameter.fs)")
    eval_string(s, "H.SamplesPerFrame = $(parameter.sample)")
    eval_string(s, "H.BitDepth = '16-bit integer' ")
    eval_string(s, "H.OutputDataType = 'double' ")
    return s
end

function (hardware::MicrophoneParameter)(select::Symbol)
    if select == :MATLAB
        return init_matlab(hardware)
    elseif select == :C++
        error("还没实现C++方法")
    elseif select == :RandomData
        error("还没实现随机数生成方法")
    else
        error("这儿没有实现其他方法")
    end
end


function getdata(s::MSession;spin = false, transform_matrix = [2, 1, 3, 4],
    ishilbert = false)
    eval_string(s, "data = H()")
    data = jarray(get_mvariable(s, :data))
    if spin
        data = view(data, :, transform_matrix)
    end
    if ishilbert
        data = hilbert(data)
    end
    return view(data, :, :)
end
