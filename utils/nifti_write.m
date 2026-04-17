function nifti_write(fname, data, hdr)
% nifti_write - 将数据写出为 NIfTI-1 单文件格式（.nii）
% 完全 standalone 实现，不依赖 SPM 或其他工具箱
%
% 输入:
%   fname - 输出文件路径（.nii）
%   data  - 数值数组，维度 [nx ny nz] 或 [nx ny nz nt]
%   hdr   - NIfTI 头结构体（由 nifti_read 返回，或手动构建）
%           必需字段: hdr.affine (4×4), hdr.pixdim (1×8 或更长)
%           可选字段: hdr.datatype（默认 16=float32）
%
% 使用示例:
%   nifti_write('output.nii', data, hdr);
%   % 从头构建:
%   hdr = nifti_default_hdr([64 64 36 200], [3 3 3 2]);
%   nifti_write('func.nii', data4d, hdr);

% -------- 强制转为 .nii --------
[fdir, fbase, fext] = fileparts(fname);
if ~strcmpi(fext, '.nii')
    fname = fullfile(fdir, [fbase '.nii']);
end

% -------- 准备维度信息 --------
siz = size(data);
while numel(siz) < 4
    siz(end+1) = 1;
end
nx = siz(1); ny = siz(2); nz = siz(3); nt = siz(4);

% -------- 确保 hdr 字段完整 --------
if ~isfield(hdr, 'datatype') || isempty(hdr.datatype)
    hdr.datatype = 16;  % float32
end
if ~isfield(hdr, 'scl_slope')
    hdr.scl_slope = 1;
end
if ~isfield(hdr, 'scl_inter')
    hdr.scl_inter = 0;
end
if ~isfield(hdr, 'descrip')
    hdr.descrip = 'Created by nifti_write';
end

% -------- 确定数据类型 --------
[writeType, bitpix, castFun] = datatype2write(hdr.datatype);

% -------- 构建 NIfTI-1 头（348字节）--------
fid = fopen(fname, 'wb', 'ieee-le');
if fid == -1
    error('[nifti_write] 无法创建文件: %s', fname);
end

% 1. sizeof_hdr (4)
fwrite(fid, int32(348), 'int32');
% 2. data_type (10) + db_name (18) = 28 bytes
fwrite(fid, zeros(1,28,'uint8'), 'uint8');
% 3. extents (4)
fwrite(fid, int32(0), 'int32');
% 4. session_error (2)
fwrite(fid, int16(0), 'int16');
% 5. regular (1)
fwrite(fid, uint8('r'), 'uint8');
% 6. dim_info (1)
fwrite(fid, uint8(0), 'uint8');

% 7. dim[8] (16)
dim = int16([4, nx, ny, nz, nt, 1, 1, 1]);
fwrite(fid, dim, 'int16');

% 8. intent_p1,p2,p3 (12)
fwrite(fid, single([0 0 0]), 'float32');
% 9. intent_code (2)
fwrite(fid, int16(0), 'int16');
% 10. datatype (2)
fwrite(fid, int16(hdr.datatype), 'int16');
% 11. bitpix (2)
fwrite(fid, int16(bitpix), 'int16');
% 12. slice_start (2)
fwrite(fid, int16(0), 'int16');

% 13. pixdim[8] (32)
if isfield(hdr, 'pixdim') && numel(hdr.pixdim) >= 4
    pd = single(hdr.pixdim(:)');
    if numel(pd) < 8
        pd(end+1:8) = 0;
    end
else
    % pixdim 从 affine 矩阵提取
    voxSz = sqrt(sum(hdr.affine(1:3,1:3).^2, 1));
    pd = single([1, voxSz(1), voxSz(2), voxSz(3), 0, 0, 0, 0]);
end
fwrite(fid, pd(1:8), 'float32');

% 14. vox_offset (4) = 352（单文件头大小 + 4字节扩展）
fwrite(fid, single(352), 'float32');
% 15. scl_slope (4)
fwrite(fid, single(hdr.scl_slope), 'float32');
% 16. scl_inter (4)
fwrite(fid, single(hdr.scl_inter), 'float32');
% 17. slice_end (2)
fwrite(fid, int16(0), 'int16');
% 18. slice_code (1)
fwrite(fid, uint8(0), 'uint8');
% 19. xyzt_units (1): mm + sec = 0x0A = 10
fwrite(fid, uint8(10), 'uint8');
% 20. cal_max, cal_min (8)
fwrite(fid, single([0 0]), 'float32');
% 21. slice_duration, toffset (8)
fwrite(fid, single([0 0]), 'float32');
% 22. glmax, glmin (8)
fwrite(fid, int32([0 0]), 'int32');

% 23. descrip (80)
desc = zeros(1,80,'uint8');
s = uint8(hdr.descrip);
n = min(numel(s), 79);
desc(1:n) = s(1:n);
fwrite(fid, desc, 'uint8');

% 24. aux_file (24)
fwrite(fid, zeros(1,24,'uint8'), 'uint8');

% 25. qform_code (2) = 0（不使用 qform）
fwrite(fid, int16(0), 'int16');
% 26. sform_code (2) = 1（使用 sform）
fwrite(fid, int16(1), 'int16');

% 27. quatern_b,c,d (12) 置零
fwrite(fid, single([0 0 0]), 'float32');
% 28. qoffset_x,y,z (12) 置零
fwrite(fid, single([0 0 0]), 'float32');

% 29. srow_x[4], srow_y[4], srow_z[4] (48)
A = hdr.affine;
fwrite(fid, single(A(1,1:4)), 'float32');
fwrite(fid, single(A(2,1:4)), 'float32');
fwrite(fid, single(A(3,1:4)), 'float32');

% 30. intent_name (16)
fwrite(fid, zeros(1,16,'uint8'), 'uint8');
% 31. magic (4): 'n+1\0'
fwrite(fid, uint8(['n','+','1',0]), 'uint8');

% -------- 4字节扩展块（必须）--------
fwrite(fid, zeros(1,4,'uint8'), 'uint8');

% -------- 写出图像数据 --------
% 如有缩放，先逆缩放
if hdr.scl_slope ~= 0 && hdr.scl_slope ~= 1
    data = (data - hdr.scl_inter) / hdr.scl_slope;
end

fwrite(fid, castFun(data(:)), writeType);
fclose(fid);

fprintf('[nifti_write] 已写出: %s  [%d %d %d %d]\n', fname, nx, ny, nz, nt);
end

% ======================================================================
% 内部辅助函数
% ======================================================================

function [typeStr, bitpix, castFun] = datatype2write(code)
switch code
    case 2
        typeStr = 'uint8';  bitpix = 8;   castFun = @uint8;
    case 4
        typeStr = 'int16';  bitpix = 16;  castFun = @int16;
    case 8
        typeStr = 'int32';  bitpix = 32;  castFun = @int32;
    case 16
        typeStr = 'float32'; bitpix = 32; castFun = @single;
    case 64
        typeStr = 'float64'; bitpix = 64; castFun = @double;
    case 256
        typeStr = 'int8';   bitpix = 8;   castFun = @int8;
    case 512
        typeStr = 'uint16'; bitpix = 16;  castFun = @uint16;
    otherwise
        typeStr = 'float32'; bitpix = 32; castFun = @single;
end
end
