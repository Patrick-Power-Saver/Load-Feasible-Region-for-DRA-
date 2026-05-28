function [NEWTON,FDFP_XB,FDFP_BX,DCPF] = idx_powerflow;

%IDX_POWERFLOW   Defines variables for power flow method.
%   [NEWTON,FDFP_XB,FDFP_BX,DCPF] = idx_powerflow;


%   by ZhaoYuan


%% define the indices
NEWTON   		= 1;	%% bus number节点编号
FDFP_XB     	= 2;	%% Pmax, maximum real power output (MW)最大有功输出
FDFP_BX 	    = 3;	%% Pmin, minimum real power output (MW)最小有功输出
DCPF        	= 4;	%% Qmax, maximum reactive power output (MVAR)最大有功输出


