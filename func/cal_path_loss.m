function [PL_dB, PL_linear, params] = cal_path_loss(distance, fc, scenario, varargin)
    % 计算路径损耗系数的统一函数
    % 输入参数:
    %   distance - 传播距离 (m), 标量或向量, > 0
    %   fc - 载波频率 (Hz), 标量, > 0
    %   scenario - 场景类型字符串 ('free_space', 'urban', 'indoor', 'suburban', 'rural')
    %   varargin - 可选参数: 'd0', 'Xsigma', 'temp', 'humidity'
    %
    % 输出参数:
    %   PL_dB - 路径损耗 (dB), 同distance维度
    %   PL_linear - 线性路径损耗系数, 同distance维度
    %   params - 参数结构体
    
    % 输入参数验证
    if nargin < 3
        error('错误: 至少需要3个输入参数 (distance, fc, scenario)');
    end
    
    % 验证distance
    if ~isnumeric(distance) || any(distance <= 0)
        error('错误: distance必须是正数');
    end
    if ~isvector(distance) && ~isscalar(distance)
        error('错误: distance必须是标量或向量');
    end
    
    % 验证fc
    if ~isnumeric(fc) || fc <= 0 || ~isscalar(fc)
        error('错误: fc必须是正标量');
    end
    
    % 验证scenario
    if ~ischar(scenario) && ~isstring(scenario)
        error('错误: scenario必须是字符串');
    end
    valid_scenarios = {'free_space', 'urban', 'indoor', 'suburban', 'rural'};
    if ~ismember(lower(scenario), valid_scenarios)
        error('错误: scenario必须是以下之一: %s', strjoin(valid_scenarios, ', '));
    end
    
    % 解析可选参数
    p = inputParser;
    addParameter(p, 'd0', 1, @(x) isnumeric(x) && x > 0);
    addParameter(p, 'Xsigma', 0, @(x) isnumeric(x) && x >= 0);
    addParameter(p, 'temp', 20, @(x) isnumeric(x) && x >= -50 && x <= 50);
    addParameter(p, 'humidity', 50, @(x) isnumeric(x) && x >= 0 && x <= 100);
    
    if nargin > 3
        parse(p, varargin{:});
    else
        parse(p);
    end
    
    d0 = p.Results.d0;
    Xsigma = p.Results.Xsigma;
    temp = p.Results.temp;
    humidity = p.Results.humidity;
    
    % 常数定义
    c = 3e8; % 光速 (m/s)
    lambda = c / fc; % 波长 (m)
    
    % 计算基本路径损耗
    scenario_lower = lower(scenario);
    
    switch scenario_lower
        case 'free_space'
            if any(distance < d0)
                warning('警告: 部分距离小于参考距离d0=%g m，可能不适用自由空间模型', d0);
            end
            PL_dB = 20*log10(4*pi*d0/lambda) + 20*log10(distance/d0);
            path_loss_exp = 2;
            
        case 'urban'
            % 3GPP Urban Macro模型
            d_break = 4 * d0 * fc / c;
            PL_dB = zeros(size(distance));
            for i = 1:length(distance)
                if distance(i) <= d_break
                    PL_dB(i) = 20*log10(4*pi*d0*fc/c) + 20*log10(distance(i)/d0);
                else
                    PL_dB(i) = 20*log10(4*pi*d0*fc/c) + 20*log10(d_break/d0) + ...
                              35*log10(distance(i)/d_break);
                end
            end
            path_loss_exp = [2, 3.5]; % 分段指数
            
        case 'indoor'
            % ITU-R P.1238室内模型
            PL_dB = 20*log10(4*pi*d0*fc/c) + 32*log10(distance/d0) + 10;
            path_loss_exp = 3.2;
            
        case 'suburban'
            % Okumura-Hata郊区模型
            d_km = distance / 1000;
            if any(d_km <= 0.1)
                warning('警告: 部分距离小于0.1km，郊区模型可能不适用');
            end
            PL_dB = 46.3 + 33.9*log10(fc/1e6) - 13.82*log10(30) + ...
                   (44.9 - 6.55*log10(30))*log10(d_km) - (3.2*(log10(11.75))^2 - 4.97);
            path_loss_exp = 3.5;
            
        case 'rural'
            % 简化农村模型
            PL_dB = 20*log10(4*pi*d0*fc/c) + 27.8*log10(distance/d0) - 20;
            path_loss_exp = 2.78;
            
        otherwise
            error('错误: 不支持的场景类型 %s', scenario);
    end
    
    % 添加阴影衰落
    if Xsigma > 0
        shadow_fading = normrnd(0, Xsigma, size(distance));
        PL_dB = PL_dB + shadow_fading;
    end
    
    % 计算线性路径损耗系数
    PL_linear = 10.^(-PL_dB/10);
    
    % 参数结构体
    params = struct();
    params.fc = fc;
    params.lambda = lambda;
    params.scenario = scenario;
    params.path_loss_exponent = path_loss_exp;
    params.reference_distance = d0;
    params.shadow_fading_std = Xsigma;
    params.temperature = temp;
    params.humidity = humidity;
    
    % 输入验证完成，返回结果
end