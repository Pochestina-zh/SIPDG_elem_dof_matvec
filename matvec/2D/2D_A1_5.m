% matrix size: 4096x4096, nnz: 64000
2D_A1_5_fid = fopen('2D_A1_5.m.dat', 'r');
2D_A1_5=fread(2D_A1_5_fid, 'double');
fclose(2D_A1_5_fid);
2D_A1_5=reshape(2D_A1_5, 3, round(length(2D_A1_5)/3));
2D_A1_5=spconvert(2D_A1_5');
