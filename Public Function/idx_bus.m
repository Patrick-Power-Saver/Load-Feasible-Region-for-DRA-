function [PQ, PV, REF, NOTE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
	VA, BASE_KV, VMAX, VMIN, PGMAX,PGMIN,QGMAX,QGMIN,PG,QG,LOLP,LOLE,LOLF,LOLD,EDNS,EENS,...
    LAM_P,LAM_Q,MU_VMAX,MU_VMIN,MU_PMAX,MU_PMIN,MU_QMAX,MU_QMIN] = idx_bus
%IDX_BUS   Defines variables for column indices to bus.
%VA, BASE_KV, VMAX, VMIN, PGMax,PGMin,QGMax,QGMin,PG,QG,LOLP,LOLE,LOLF,LOLD,EDNS,EENS,...
%    LAM_P,LAM_Q,MU_VMAX,MU_VMIN,MU_PMAX,MU_PMIN,MU_QMAX,MU_QMIN] = idx_bus

%   by ZhaYuan


%% define bus types节点类型
PQ		= 1;
PV		= 2;
REF		= 3;
NOTE    = 4; 

%% define the indices
BUS_I		= 1;	%% bus number (1 to 29997)节点编号
BUS_TYPE	= 2;	%% bus type (1 - PQ bus, 2 - PV bus, 3 - reference bus, 4 - isolated bus)节点类型
PD			= 3;	%% Pd, real power demand (MW)有功需求
QD			= 4;	%% Qd, reactive power demand (MVAR)无功需求
GS			= 5;	%% Gs, shunt conductance (MW at V = 1.0 p.u.)并联电导
BS			= 6;	%% Bs, shunt susceptance (MVAR at V = 1.0 p.u.)并联电纳
BUS_AREA	= 7;	%% area number, 1-100区域编号
VM			= 8;	%% Vm, voltage magnitude (p.u.)电压幅值
VA			= 9;	%% Va, voltage angle (degrees)电压相角（度）
BASE_KV		= 10;	%% baseKV, base voltage (kV)参考电压
VMAX		= 11;	%% maxVm, maximum voltage magnitude (p.u.)		(not in PTI format)最大电压幅值
VMIN		= 12;	%% minVm, minimum voltage magnitude (p.u.)		(not in PTI format)最小电压幅值

PGMAX       = 13;   %% maximum active power wthich can be ejected into the net by the generator installed on this bus最大发出有功
PGMIN       = 14;   %% PGMin, minimum active power wthich can be ejected into the net by the generator installed on this bus最小发出有功
QGMAX       = 15;   %% maximum reactive power wthich can be ejected into the net by the generator installed on this bus最大发出无功
QGMIN       = 16;   %% QGMin, minimum reactive power wthich can be ejected into the net by the generator installed on this bus最小发出有功
PG          = 17;   %% PG, actual active power wthich is ejected into the net by the generator installed on this bus实际注入有功
QG          = 18;   %% actual reactive power wthich is ejected into the net by the generator installed on this bus实际注入无功
LOLP        = 19;   %% LOLP, Loss of Load Probability缺电概率
LOLE        = 20;   %% LOLE, Loss of Load Expection缺电时间期望
LOLF        = 21;   %% Loss of Load Frequency缺电频率
LOLD        = 22;   %% LOLD, Loss of Load Duration缺电持续时间
EDNS        = 23;   %% Expected Demand Not Supplied期望缺供电力
EENS        = 24;   %% Expected Energy Not Supplied期望缺共电量
LAM_P       = 25;   %% Lagrange multiplier on real power mismatch拉格朗日乘子 
LAM_Q       = 26;   %% Lagrange multiplier on reactive power mismatch
MU_VMAX     = 27;   %% Kuhn-Tucker multiplier on upper voltage limit
MU_VMIN     = 28;   %% Kuhn-Tucker multiplier on lower voltage limit
MU_PMAX     = 29;   %% Kuhn-Tucker multiplier on upper Pg limit
MU_PMIN     = 30;   %% Kuhn-Tucker multiplier on lower Pg limit
MU_QMAX     = 31;   %% Kuhn-Tucker multiplier on upper Qg limit
MU_QMIN     = 32;   %% Kuhn-Tucker multiplier on lower Qg limit
return;
