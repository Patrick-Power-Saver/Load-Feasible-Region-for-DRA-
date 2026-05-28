function [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE, TAP, SHIFT_ANGLE, BR_STATUS, BR_LAMDA,...
        BR_MTTF,TCSC_X, BR_PFAILURE, PF, QF, PT, QT, MU_SF, MU_ST, KT_MAXANGLE_TPST,KT_MINANGLE_TPST,KT_MAX_X_TCSC,KT_MIN_X_TCSC] = idx_brch
%IDX_BRCH   Defines variables for column indices to branch.
%[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE, ...
%TAP, SHIFT_ANGLE, BR_STATUS, BR_LAMDA, BR_MTTF, BR_PFAILURE, PF, QF, PT, QT, MU_SF, MU_ST, KT_ANGLE] = idx_brch

%   by ZhaoYuan


%% define the indices
F_BUS		= 1;	%% f, from bus number起始节点编号
T_BUS		= 2;	%% t, to bus number终止节点编号
BR_R		= 3;	%% r, resistance (p.u.)电阻
BR_X		= 4;	%% x, reactance (p.u.)电抗
BR_B		= 5;	%% b, total line charging susceptance (p.u.)充电电纳
RATE		= 6;	%% rateA, MVA rating A (long term rating)
TAP			= 7;	%% ratio, transformer off nominal turns ratio变压器非标准变比
SHIFT_ANGLE	= 8;	%% angle, transformer phase shift angle变压器移相角
BR_STATUS	= 9;	%% initial branch status, 1 - in service, 0 - out of service支路状态1为运行（0为故障）
BR_LAMDA    = 10;   %% Failure rate of transmission line变压器支路故障率
BR_MTTF     = 11;   %% Repair duration修复时间      
TCSC_X      = 12;   %% TCSC Xc Value可控串联补偿器

BR_PFAILURE = 13;   %% Failure Probability故障概率
%% included in power flow solution and optimal load shedding, not
%% necessarily in input在潮流计算和最优负荷削减中需要的数据，输入中不需要
PF			= 14;	%% real power injected at "from" bus end (MW)		(not in PTI format)始节点注入有功
QF			= 15;	%% reactive power injected at "from" bus end (MVAR)	(not in PTI format)终节点注入有功
PT			= 16;	%% real power injected at "to" bus end (MW)			(not in PTI format)始节点注入无功
QT			= 17;	%% reactive power injected at "to" bus end (MVAR)	(not in PTI format)终节点注入无功

%% included in opf solution, not necessarily in
%% input包括在最优负荷削减中，在输入数据中没有必要输入
%% assume objective function has units, u
MU_SF		= 18;	%% Kuhn-Tucker multiplier on MVA limit at "from" bus (u/MVA)
MU_ST		= 19;	%% Kuhn-Tucker multiplier on MVA limit at "to" bus (u/MVA)
KT_MAXANGLE_TPST = 20;   %% Kuhn-Tucker multiplier on shifting transformer angle limit 
KT_MINANGLE_TPST = 21;   %% Kuhn-Tucker multiplier on shifting transformer angle limit
KT_MAX_X_TCSC = 22;   %%Kuhn-Tucker multiplier on TSCS Xc upper  limit
KT_MIN_X_TCSC = 23;   %%Kuhn-Tucker multiplier on TCSC Xc lower  limit

return;





