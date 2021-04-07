
fid=fopen('data.txt');
for k=1:9
fgetl(fid);
end
C = textscan(fid,'%s = %s');
fclose(fid);



fid=fopen('../sim/circuit1.txt', 'wt');

fprintf(fid, 'Vs 1 0 DC %s\n', C{2}{8});
fprintf(fid, 'R1 1 2 %sK\n', C{2}{1});
fprintf(fid, 'R2 2 3 %sK\n', C{2}{2});
fprintf(fid, 'R3 5 2 %sK\n', C{2}{3});
fprintf(fid, 'R4 5 0 %sK\n', C{2}{4});
fprintf(fid, 'R5 6 5 %sK\n', C{2}{5});
fprintf(fid, 'R6 0 7 %sK\n', C{2}{6});
fprintf(fid, 'R7 9 8 %sK\n', C{2}{7});
fprintf(fid, 'GIb 6 3 (2,5) %sm\n', C{2}{10});
fprintf(fid, 'VIc 7 9 0V\n');
fprintf(fid, 'HVd 5 8 VIc %sK\n', C{2}{11});
fprintf(fid, 'Cap 6 8 %sUF\n', C{2}{9});

fclose(fid);


Vx = V6 - V8;

fid=fopen('data.txt');
for k=1:9
fgetl(fid);
end
C = textscan(fid,'%s = %s');
fclose(fid);

fid=fopen('../sim/circuit2.txt','wt');

fprintf(fid, 'Vs 1 0 0V\n');
fprintf(fid, 'R1 1 2 %sK\n', C{2}{1});
fprintf(fid, 'R2 2 3 %sK\n', C{2}{2});
fprintf(fid, 'R3 5 2 %sK\n', C{2}{3});
fprintf(fid, 'R4 5 0 %sK\n', C{2}{4});
fprintf(fid, 'R5 6 5 %sK\n', C{2}{5});
fprintf(fid, 'R6 0 7 %sK\n', C{2}{6});
fprintf(fid, 'R7 8 9 %sK\n', C{2}{7});
fprintf(fid, 'GIb 6 3 (2,5) %sm\n', C{2}{10});
fprintf(fid, 'VIc 7 9 0V\n');
fprintf(fid, 'HVd 5 8 VIc %sK\n', C{2}{11});
fprintf(fid, 'Vx 6 8 %.15f', Vx);

fclose(fid);
