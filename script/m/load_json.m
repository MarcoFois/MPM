function stct = load_json (filename)
  fid = fopen (filename, 'r');
  s = fread (fid, inf, 'char');
  fclose (fid);
  stct = jsondecode (char(s).');
endfunction

