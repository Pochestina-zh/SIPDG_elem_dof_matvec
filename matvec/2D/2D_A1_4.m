% matrix size: 1024x1024, nnz: 15616
2D_A1_4_fid = fopen('2D_A1_4.m.dat', 'r');
2D_A1_4=fread(2D_A1_4_fid, 'double');
fclose(2D_A1_4_fid);
2D_A1_4=reshape(2D_A1_4, 3, round(length(2D_A1_4)/3));
2D_A1_4=spconvert(2D_A1_4');
