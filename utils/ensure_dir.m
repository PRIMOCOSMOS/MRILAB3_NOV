function ensure_dir(dirPath)
% ensure_dir - 确保目录存在，若不存在则递归创建
%
% 输入:
%   dirPath - 目录路径字符串
%
% 使用示例:
%   ensure_dir('D:\MRI_PRO\MRILAB3\output\Sub_01');

if ~exist(dirPath, 'dir')
    [status, msg] = mkdir(dirPath);
    if status
        fprintf('[ensure_dir] 已创建目录: %s\n', dirPath);
    else
        error('[ensure_dir] 无法创建目录 %s: %s', dirPath, msg);
    end
end
end
