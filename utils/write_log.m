function write_log(logFile, msg)
% write_log - 将带时间戳的消息写入日志文件并同步打印到命令窗口
%
% 输入:
%   logFile - 日志文件完整路径（字符串）
%   msg     - 消息字符串
%
% 使用示例:
%   write_log('/data/logs/pipeline.log', 'Step01: DICOM 转换完成');

timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');
line = sprintf('[%s] %s\n', timestamp, msg);

% 打印到命令窗口
fprintf('%s', line);

% 追加写入日志文件
fid = fopen(logFile, 'a');
if fid ~= -1
    fprintf(fid, '%s', line);
    fclose(fid);
else
    warning('[write_log] 无法打开日志文件: %s', logFile);
end
end
