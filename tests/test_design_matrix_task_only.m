function tests = test_design_matrix_task_only
% test_design_matrix_task_only - 回归测试：单任务范式应使用隐式 Rest 基线
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
rootDir = fileparts(fileparts(mfilename('fullpath')));
addpath(rootDir);
addpath(fullfile(rootDir, 'stats'));

testCase.TestData.rootDir = rootDir;
end

function test_config_uses_single_task_condition(testCase)
cfg = config_sub01();

verifyEqual(testCase, numel(cfg.cond.names), 1);
verifyEqual(testCase, cfg.cond.names{1}, 'Task');
verifyEqual(testCase, numel(cfg.cond.onsets), 1);
verifyEqual(testCase, numel(cfg.cond.durations), 1);

verifyEqual(testCase, numel(cfg.tcons.weight), 1);
verifyEqual(testCase, cfg.tcons.weight(1), 1);
end

function test_design_matrix_has_one_task_column(testCase)
cfg = config_sub01();
nScans = 150;

[X, names] = build_design_matrix(cfg, nScans);

verifyEqual(testCase, size(X, 1), nScans);
verifyEqual(testCase, names{1}, 'Task');
verifyFalse(testCase, any(strcmp(names, 'Rest')));
verifyGreaterThan(testCase, std(X(:,1)), 1e-6);

nDrift = fix(2 * nScans * cfg.TR / cfg.hpf + 1) - 1;
useExplicitDrift = true;
if isfield(cfg, 'glm') && isfield(cfg.glm, 'explicitDriftRegressors')
	useExplicitDrift = logical(cfg.glm.explicitDriftRegressors);
end
if useExplicitDrift
	expectedCols = 1 + 1 + max(0, nDrift);
else
	expectedCols = 1 + 1;
end
verifyEqual(testCase, size(X, 2), expectedCols);
end
