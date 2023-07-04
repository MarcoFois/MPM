function save_json (stct, filename, prettyprint = false)
  s = jsonencode (stct, 'prettyprint', prettyprint);
  fid = fopen (filename, 'w');
  fprintf (fid, s);
  fclose (fid);
endfunction
