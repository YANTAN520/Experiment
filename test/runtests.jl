using Experiment
using Test

@testset "Experiment.jl" begin
    ## 测试读取数据是否出错
    M = MicrophoneParameter("线路输入 (LIS_MIIS_8X8_192K)", 44100, 30, 4);
    hardware = M(:MATLAB);
    data = getdata(hardware, ishilbert = true);

    ## 测试过备原子基生成, 全部ok，但是二维会有181*181个基
    ##窄带情况都无法满足，宽带情况更加无法时间要求
    l1 = linear1(0.17, 4);
    l2 = linear2(0.17, 0.1, 4, 4);
    C = circle(4, 0.17);

    Al1 = A_matrix((0.0, 90.0), [0.17], 0.5*pi/180, l1, 1);
    Al2 = A_matrix((0.0, 90.0, 0.0, 90.0), [0.17], [0.5*pi/180; 0.5*pi/180], l2, 2);
    Al3 = A_matrix((0.0, 90.0, 0.0, 90.0), [0.17], [0.5*pi/180; 0.5*pi/180], C, 2);

    A1 = Al1();
    A2 = Al2();
    A3 = Al3();
    ##函数测试全部ok
    x = randn(1000,10);
    xx = add_gaussian_noise(x, 20.0, dims = 1);
    xxx = hilbert(x, dims = 1);

    find_fft_peak(x[:,1], 44100, 0.01);
    find_fft_peak(x[:,2], 44100, 3);

    ##测试solver函数
    P1 = Position(data', A1, 1, :svd_l1!);
    # P2 = Position(data', A2, 2, :svd_l1!);
    P3 = Position(data', A3, 2, :svd_l1!);

    solver(P1,true);
    solver(P3, true);

end

@testset "plot position" begin
    
end

