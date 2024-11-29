% matrix size: 262144x262144, nnz: 4.18202e+06
2D_A1_8_fid = fopen('2D_A1_8.m.dat', 'r');
2D_A1_8=fread(2D_A1_8_fid, 'double');
fclose(2D_A1_8_fid);
2D_A1_8=reshape(2D_A1_8, 3, round(length(2D_A1_8)/3));
2D_A1_8=spconvert(2D_A1_8');
