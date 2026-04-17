function [data, hdr] = nifti_read(fname)
% nifti_read - 读取 NIfTI-1 格式文件（.nii 单文件或 .hdr/.img 双文件）
% 完全 standalone 实现，不依赖 SPM 或其他工具箱
%
% 输入:
%   fname - NIfTI 文件路径（.nii 或 .hdr，支持 .nii.gz 需要系统有 gunzip）
%
% 输出:
%   data  - 数值数组，维度为 [nx ny nz nt]（nt=1 时为3D）
%   hdr   - NIfTI 头结构体（见下方字段说明）
%
% 支持的数据类型: uint8, int16, int32, float32, float64, int8, uint16, uint32
%
% 使用示例:
%   [data, hdr] = nifti_read('func.nii');
%   [data, hdr] = nifti_read('anat.hdr');

% -------- 处理 .gz 压缩 --------
[fdir, fbase, fext] = fileparts(fname);
if strcmpi(fext, '.gz')
    tmpFile = fullfile(tempdir, fbase);
    gunzip(fname, tempdir);
    fname = tmpFile;
    [~, fbase, fext] = fileparts(fname);
end

% -------- 确定头文件与数据文件 --------
if strcmpi(fext, '.hdr')
    hdrFile = fname;
    imgFile = fullfile(fdir, [fbase '.img']);
    isSingleFile = false;
elseif strcmpi(fext, '.nii')
    hdrFile = fname;
    imgFile = fname;
    isSingleFile = true;
else
    error('[nifti_read] 不支持的文件扩展名: %s', fext);
end

% -------- 打开并读取头信息 --------
fid = fopen(hdrFile, 'rb', 'ieee-le');
if fid == -1
    error('[nifti_read] 无法打开文件: %s', hdrFile);
end

% 先读 sizeof_hdr 判断字节序
raw4 = fread(fid, 4, 'uint8');
sizeof_hdr = typecast(uint8(raw4), 'int32');
if sizeof_hdr ~= 348
    % 尝试大端
    fclose(fid);
    fid = fopen(hdrFile, 'rb', 'ieee-be');
    raw4 = fread(fid, 4, 'uint8');
    sizeof_hdr = typecast(uint8(raw4), 'int32');
    if sizeof_hdr ~= 348
        error('[nifti_read] 无效的 NIfTI 文件（sizeof_hdr != 348）: %s', hdrFile);
    end
    byteOrder = 'ieee-be';
else
    byteOrder = 'ieee-le';
    fclose(fid);
    fid = fopen(hdrFile, 'rb', byteOrder);
end

% -------- 解析 NIfTI-1 头（共 348 字节）--------
hdr.sizeof_hdr    = fread(fid, 1, 'int32');        % 4
hdr.data_type     = deblank(fread(fid, 10, '*char')');  % 10
hdr.db_name       = deblank(fread(fid, 18, '*char')');  % 18
hdr.extents       = fread(fid, 1, 'int32');        % 4
hdr.session_error = fread(fid, 1, 'int16');        % 2
hdr.regular       = fread(fid, 1, '*char');        % 1
hdr.dim_info      = fread(fid, 1, 'uint8');        % 1

hdr.dim           = fread(fid, 8, 'int16');        % 16
hdr.intent_p1     = fread(fid, 1, 'float32');      % 4
hdr.intent_p2     = fread(fid, 1, 'float32');      % 4
hdr.intent_p3     = fread(fid, 1, 'float32');      % 4
hdr.intent_code   = fread(fid, 1, 'int16');        % 2
hdr.datatype      = fread(fid, 1, 'int16');        % 2
hdr.bitpix        = fread(fid, 1, 'int16');        % 2
hdr.slice_start   = fread(fid, 1, 'int16');        % 2
hdr.pixdim        = fread(fid, 8, 'float32');      % 32
hdr.vox_offset    = fread(fid, 1, 'float32');      % 4
hdr.scl_slope     = fread(fid, 1, 'float32');      % 4
hdr.scl_inter     = fread(fid, 1, 'float32');      % 4
hdr.slice_end     = fread(fid, 1, 'int16');        % 2
hdr.slice_code    = fread(fid, 1, 'uint8');        % 1
hdr.xyzt_units    = fread(fid, 1, 'uint8');        % 1
hdr.cal_max       = fread(fid, 1, 'float32');      % 4
hdr.cal_min       = fread(fid, 1, 'float32');      % 4
hdr.slice_duration= fread(fid, 1, 'float32');      % 4
hdr.toffset       = fread(fid, 1, 'float32');      % 4
hdr.glmax         = fread(fid, 1, 'int32');        % 4
hdr.glmin         = fread(fid, 1, 'int32');        % 4
hdr.descrip       = deblank(fread(fid, 80, '*char')');  % 80
hdr.aux_file      = deblank(fread(fid, 24, '*char')');  % 24
hdr.qform_code    = fread(fid, 1, 'int16');        % 2
hdr.sform_code    = fread(fid, 1, 'int16');        % 2
hdr.quatern_b     = fread(fid, 1, 'float32');      % 4
hdr.quatern_c     = fread(fid, 1, 'float32');      % 4
hdr.quatern_d     = fread(fid, 1, 'float32');      % 4
hdr.qoffset_x     = fread(fid, 1, 'float32');      % 4
hdr.qoffset_y     = fread(fid, 1, 'float32');      % 4
hdr.qoffset_z     = fread(fid, 1, 'float32');      % 4
hdr.srow_x        = fread(fid, 4, 'float32');      % 16
hdr.srow_y        = fread(fid, 4, 'float32');      % 16
hdr.srow_z        = fread(fid, 4, 'float32');      % 16
hdr.intent_name   = deblank(fread(fid, 16, '*char')');  % 16
hdr.magic         = deblank(fread(fid, 4, '*char')');   % 4
fclose(fid);

% -------- 构建仿射矩阵 --------
% 优先使用 sform（sform_code > 0）
if hdr.sform_code > 0
    hdr.affine = [hdr.srow_x(:)'; hdr.srow_y(:)'; hdr.srow_z(:)'; 0 0 0 1];
elseif hdr.qform_code > 0
    hdr.affine = qform2affine(hdr);
else
    % 使用 pixdim 构建对角仿射
    hdr.affine = diag([hdr.pixdim(2:4); 1]);
end

% -------- 提取维度 --------
ndim = hdr.dim(1);
hdr.nx = hdr.dim(2);
hdr.ny = hdr.dim(3);
hdr.nz = max(hdr.dim(4), 1);
hdr.nt = max(hdr.dim(5), 1);
hdr.byteOrder = byteOrder;

% -------- 读取图像数据 --------
% 确定跳过的字节数
if isSingleFile
    offset = max(hdr.vox_offset, 352);
else
    offset = 0;
end

fid2 = fopen(imgFile, 'rb', byteOrder);
if fid2 == -1
    error('[nifti_read] 无法打开数据文件: %s', imgFile);
end
fseek(fid2, offset, 'bof');

% 确定读取类型
[readType, nBytes] = datatype2str(hdr.datatype);
nVox = hdr.nx * hdr.ny * hdr.nz * hdr.nt;
raw = fread(fid2, nVox, readType);
fclose(fid2);

if numel(raw) < nVox
    warning('[nifti_read] 实际读取 %d 体素，期望 %d，自动补零', numel(raw), nVox);
    raw(end+1:nVox) = 0;
end

% 重塑为4D矩阵
data = reshape(raw, hdr.nx, hdr.ny, hdr.nz, hdr.nt);
data = double(data);

% 应用缩放因子
if hdr.scl_slope ~= 0 && ~isnan(hdr.scl_slope)
    data = data * hdr.scl_slope + hdr.scl_inter;
end

end % function nifti_read

% ======================================================================
% 内部辅助函数
% ======================================================================

function [typeStr, nBytes] = datatype2str(code)
% 将 NIfTI datatype 代码映射到 MATLAB fread 类型字符串
switch code
    case 2,   typeStr = 'uint8';   nBytes = 1;
    case 4,   typeStr = 'int16';   nBytes = 2;
    case 8,   typeStr = 'int32';   nBytes = 4;
    case 16,  typeStr = 'float32'; nBytes = 4;
    case 32,  typeStr = 'float32'; nBytes = 4;   % complex64（取实部）
    case 64,  typeStr = 'float64'; nBytes = 8;
    case 256, typeStr = 'int8';    nBytes = 1;
    case 512, typeStr = 'uint16';  nBytes = 2;
    case 768, typeStr = 'uint32';  nBytes = 4;
    otherwise
        warning('[nifti_read] 未知 datatype %d，默认使用 float32', code);
        typeStr = 'float32'; nBytes = 4;
end
end

function affine = qform2affine(hdr)
% 从四元数（quaternion）和偏移量构建 4×4 仿射矩阵
% 参考 NIfTI-1 标准公式
b = hdr.quatern_b;
c = hdr.quatern_c;
d = hdr.quatern_d;
a = sqrt(max(1 - b^2 - c^2 - d^2, 0));

R = [a^2+b^2-c^2-d^2,  2*(b*c-a*d),      2*(b*d+a*c);
     2*(b*c+a*d),      a^2+c^2-b^2-d^2,  2*(c*d-a*b);
     2*(b*d-a*c),      2*(c*d+a*b),      a^2+d^2-b^2-c^2];

dx = hdr.pixdim(2);
dy = hdr.pixdim(3);
dz = hdr.pixdim(4) * sign(hdr.pixdim(1));  % qfac

T = [hdr.qoffset_x; hdr.qoffset_y; hdr.qoffset_z];

affine = [R .* [dx dy dz], T; 0 0 0 1];
end
