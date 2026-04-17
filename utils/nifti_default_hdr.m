function hdr = nifti_default_hdr(dims, voxSize)
% nifti_default_hdr - 构建最小化有效 NIfTI 头结构体
%
% 输入:
%   dims    - 图像维度 [nx ny nz] 或 [nx ny nz nt]
%   voxSize - 体素尺寸 [dx dy dz] 或 [dx dy dz dt]（mm 和 秒）
%
% 输出:
%   hdr - NIfTI 头结构体，可直接传给 nifti_write
%
% 使用示例:
%   hdr = nifti_default_hdr([64 64 36 200], [3 3 3 2]);

if nargin < 2
    voxSize = ones(1, numel(dims));
end

dims    = dims(:)';
voxSize = voxSize(:)';

while numel(dims) < 4,    dims(end+1)    = 1; end
while numel(voxSize) < 4, voxSize(end+1) = 0; end

nx = dims(1); ny = dims(2); nz = dims(3); nt = dims(4);
dx = voxSize(1); dy = voxSize(2); dz = voxSize(3); dt = voxSize(4);

hdr.sizeof_hdr    = 348;
hdr.data_type     = '';
hdr.db_name       = '';
hdr.extents       = 0;
hdr.session_error = 0;
hdr.regular       = 'r';
hdr.dim_info      = 0;
hdr.dim           = int16([4, nx, ny, nz, nt, 1, 1, 1]);
hdr.intent_p1     = 0;
hdr.intent_p2     = 0;
hdr.intent_p3     = 0;
hdr.intent_code   = 0;
hdr.datatype      = 16;    % float32
hdr.bitpix        = 32;
hdr.slice_start   = 0;
hdr.pixdim        = single([1, dx, dy, dz, dt, 0, 0, 0]);
hdr.vox_offset    = 352;
hdr.scl_slope     = 1;
hdr.scl_inter     = 0;
hdr.slice_end     = 0;
hdr.slice_code    = 0;
hdr.xyzt_units    = 10;    % mm + sec
hdr.cal_max       = 0;
hdr.cal_min       = 0;
hdr.slice_duration= 0;
hdr.toffset       = 0;
hdr.glmax         = 0;
hdr.glmin         = 0;
hdr.descrip       = 'nifti_default_hdr';
hdr.aux_file      = '';
hdr.qform_code    = 0;
hdr.sform_code    = 1;
hdr.quatern_b     = 0;
hdr.quatern_c     = 0;
hdr.quatern_d     = 0;
hdr.qoffset_x     = 0;
hdr.qoffset_y     = 0;
hdr.qoffset_z     = 0;

% 默认仿射：以图像中心为原点
hdr.affine = [dx 0  0  -dx*nx/2;
              0  dy 0  -dy*ny/2;
              0  0  dz -dz*nz/2;
              0  0  0   1      ];

hdr.srow_x = hdr.affine(1,:);
hdr.srow_y = hdr.affine(2,:);
hdr.srow_z = hdr.affine(3,:);
hdr.intent_name = '';
hdr.magic       = 'n+1';
hdr.nx = nx; hdr.ny = ny; hdr.nz = nz; hdr.nt = nt;
hdr.byteOrder = 'ieee-le';
end
