function write_case(outfile)
fprintf(outfile,'FORMAT\n\n');
fprintf(outfile,'type:\tensight gold\n\n');
fprintf(outfile,'GEOMETRY\n\n');
fprintf(outfile,'model:\tsolid.geo\n\n');
fprintf(outfile,'VARIABLE\n\n');
fprintf(outfile,'tensor symm per element:\tsigma\tsolid.sig');
end