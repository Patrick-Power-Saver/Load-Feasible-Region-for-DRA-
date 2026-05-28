%% R_main_reliability_evaluation_for Dist.m
% 基于交流负荷可行域（AC-LFR）的配电网可靠性评估主脚本
%
% ===================================================================
%  相较于 R_main_reliability_evaluation.m（DC版）的改动
% ===================================================================
%  改动1（~第42行）：建库函数
%    DC版：LFR_lib = build_LFR_library_v2(CASENAME, opts)
%    AC版：LFR_lib = build_LFR_library_v2_AC(CASENAME, opts)
%
%  改动2（~第72行）：lib_filename 区分 AC/DC
%    DC版：lib_filename = sprintf('LFR_lib_%s_%s.mat', ...)
%    AC版：lib_filename = sprintf('AC_LFR_lib_%s_%s.mat', ...)
%
%  其余所有代码（MC评估、收敛判断、结果汇总、绘图）完全不变。
% ===================================================================

clear; clc; close all;
fprintf('========================================\n');
fprintf('  基于交流负荷可行域的配电网可靠性评估\n');
fprintf('========================================\n\n');

%% ===================================================================
%  用户配置区
% ====================================================================
CASENAME = 'dist_opf_case85';%dist_opf_case85 dist_opf_case137 dist_opf_case417 dist_opf_case1080 dist_opf_case33_ieee RTS96 RTS79_1 trans_opf_case_rbts
%convert_data_to_m('85-Node System Data.xlsx', 'dist_opf_case85.m', 'dist_opf_case85');

opts.rebuild_offline_lib = true;
opts.max_samples         = 2e6;
opts.conv_eps            = 0.001;
opts.conv_check          = 2000;
opts.verbose             = true;

FIXED_SEED = 2026;

[baseMVA, bus, ~,~,~,~,~,~,~,~,~,~,~] = feval(CASENAME);
L0 = bus(:, 3) / baseMVA;

strategies = {'Original_Theta', 'Improved_PBC_Theta'};
results = struct();

fprintf('==================================================\n');
fprintf('  AC-LFR 可靠性评估\n');
fprintf('  算例: %s | 重新建库: %d | 固定种子: %d\n', ...
        CASENAME, opts.rebuild_offline_lib, FIXED_SEED);
fprintf('==================================================\n');

for i = 1:length(strategies)
    mode_name = strategies{i};
    fprintf('\n>>> 当前评估策略: [%s]\n', mode_name);

    if strcmp(mode_name, 'Original_Theta')
        opts.use_tight_theta = false;
    else
        opts.use_tight_theta = true;
    end

    % ★ 改动2：文件名加 AC_ 前缀以区分 DC 版库文件 ★
    lib_filename = sprintf('AC_LFR_lib_%s_%s.mat', CASENAME, mode_name);
    
    if opts.rebuild_offline_lib || ~exist(lib_filename, 'file')
        fprintf('  [离线阶段] 正在生成交流 LFR 库...\n');
        t_build = tic;
        % ★ 改动1：调用交流版建库函数 ★
        LFR_lib = build_LFR_library_v2_Dist(CASENAME, opts);
        time_build = toc(t_build);
        save(lib_filename, 'LFR_lib', 'time_build');
        fprintf('  [离线阶段] 建库完成 (耗时: %.2fs)\n', time_build);
    else
        fprintf('  [离线阶段] 加载已有 AC-LFR 库: %s\n', lib_filename);
        load(lib_filename, 'LFR_lib', 'time_build');
        fprintf('  [离线阶段] 加载完成 (原建库耗时: %.2fs)\n', time_build);
    end

    %% 在线蒙特卡洛评估
    fprintf('  [在线阶段] 开始蒙特卡洛抽样...\n');
    rng(FIXED_SEED);

    t_mc = tic;
    [EENS, LOLP, details] = mc_reliability_with_LFR_for_Dist_v3( ...
        CASENAME, LFR_lib, L0, opts);
    time_mc = toc(t_mc);

    fprintf('  [在线阶段] 完成 (耗时: %.2fs，迭代: %d次)\n', ...
            time_mc, details.n_sample);

    T_period = 8760;
    results.(mode_name).LOLP       = LOLP;
    results.(mode_name).EENS       = EENS * baseMVA * T_period;
    results.(mode_name).details    = details;
    results.(mode_name).time_build = time_build;
    results.(mode_name).time_mc    = time_mc;

    plot_reliability_summary(details, results.(mode_name).EENS, LOLP);
end

%% 结果对比
fprintf('\n==================================================\n');
fprintf('              AC-LFR 对比测试结果汇总\n');
fprintf('==================================================\n');
fprintf('%-20s | %-10s | %-12s | %-10s | %-10s\n', ...
    '策略', 'LOLP', 'EENS(MWh/yr)', '建库(s)', '抽样(s)');
fprintf('------------------------------------------------------------------\n');
for i = 1:length(strategies)
    m = strategies{i};
    fprintf('%-20s | %-10.6f | %-12.1f | %-10.2f | %-10.2f\n', ...
        m, results.(m).LOLP, results.(m).EENS, ...
        results.(m).time_build, results.(m).time_mc);
end

%% 绘图函数（与DC版完全相同）
function plot_reliability_summary(details, EENS_MWh, LOLP_val)
figure('Name','AC-LFR 可靠性评估结果','Position',[100 100 920 420]);

subplot(1,2,1);
times  = [details.t_adequacy, details.t_loadshed, ...
          details.t_correct,  details.t_cache];
labels = {'充裕性判断','切负荷计算','LFR修正','缓存查询'};
valid  = times > 1e-6;
if any(valid)
    pie(times(valid) + eps, labels(valid));
    title('计算时间分解（在线阶段）','FontSize',12);
else
    text(0.5,0.5,'无计时数据','HorizontalAlignment','center'); axis off;
end

subplot(1,2,2);
n_LFR_new = max(0, details.n_LFR_used - details.n_cache_hit);
counts    = [details.n_adequate, details.n_cache_hit, n_LFR_new];
labels2   = {'充裕（跳过）','缓存命中','LFR新建'};
bar(counts);
set(gca,'XTickLabel',labels2,'XTickLabelRotation',12,'FontSize',9);
ylabel('次数'); grid on;
title(sprintf('状态分布（共%d次）',details.n_sample),'FontSize',12);
sgtitle(sprintf('AC-LFR: LOLP=%.4f，EENS=%.1f MWh/年', LOLP_val, EENS_MWh), ...
        'FontSize',13,'FontWeight','bold');
end
