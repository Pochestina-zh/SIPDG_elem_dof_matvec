% matrix size: 16384x16384, nnz: 259072
2D_A1_6_fid = fopen('2D_A1_6.m.dat', 'r');
2D_A1_6=fread(2D_A1_6_fid, 'double');
fclose(2D_A1_6_fid);
2D_A1_6=reshape(2D_A1_6, 3, round(length(2D_A1_6)/3));
2D_A1_6=spconvert(2D_A1_6');
