% matrix size: 65536x65536, nnz: 1.04243e+06
2D_A1_7_fid = fopen('2D_A1_7.m.dat', 'r');
2D_A1_7=fread(2D_A1_7_fid, 'double');
fclose(2D_A1_7_fid);
2D_A1_7=reshape(2D_A1_7, 3, round(length(2D_A1_7)/3));
2D_A1_7=spconvert(2D_A1_7');
