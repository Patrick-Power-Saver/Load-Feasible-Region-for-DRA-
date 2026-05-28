function [GEN_BUS, GEN_PMAX, GEN_PMIN, GEN_QMAX, GEN_QMIN, GEN_STATUS, ...
	GEN_LAMDA,GEN_MTTF,GEN_PG,GEN_QG,GEN_PFAILURE] = idx_gen
%IDX_GEN   Defines variables for column indices to gen.
%   [GEN_BUS, GEN_PMAX, GEN_PMIN, GEN_QMAX, GEN_QMIN, GEN_STATUS, ...
%	GEN_LAMDA, GEN_PFAILURE] = idx_gen

%   MATPOWER Version 2.0
%   by Ray Zimmerman, PSERC Cornell    9/19/97
%   Copyright (c) 1996, 1997 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/ for more info.

%% define the indices
GEN_BUS 		= 1;	%% bus number节点编号
GEN_PMAX    	= 2;	%% Pmax, maximum real power output (MW)最大有功输出
GEN_PMIN	    = 3;	%% Pmin, minimum real power output (MW)最小有功输出
GEN_QMAX    	= 4;	%% Qmax, maximum reactive power output (MVAR)最大无功输出
GEN_QMIN    	= 5;	%% Qmin, minimum reactive power output (MVAR)最小无功输出
GEN_STATUS	    = 6;	%% status, 1 - machine in service, 0 - machine out of service发点击状态1为运行状态（0为故障状态）
GEN_LAMDA       = 7;    %% Failure Rate故障率     
GEN_MTTF        = 8;    %% Repair duration修复时间
GEN_PG          = 9;    %% original real power output (MW)初始有功输出
GEN_QG          = 10;   %% original reactive power output (MVAR)初始无功输出
GEN_PFAILURE    = 11;    %% PFailure, Probability of Failure故障概率



