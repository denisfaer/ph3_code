program ph3_curves;

const
  cond = 6;         // number of experimental groups INTEGER-ONLY
  k = 15;           // max stained events per gonad INTEGER-ONLY
  n = 100;          // max samples per condition INTEGER-ONLY
  rmin = 0;         // min input value
  rmax = 100;       // max input value
  fine = 1;         // sampling step (difference between adjacent X values in generated curves)
  window = 5;       // density kernel window side
  entries = 100;    // data reseving constant, please re-calculate manually: Ent = (rmax - rmin) / fine [this gives you minimal amout of values to reserve, input next highest integer]
  pairs = 105;      // data reseving constant, please re-calculate manually: Pr = k * (k - 1) / 2
  randpair = FALSE; // generate random data and compute pairs

var
  fin, fout: text;
  y, j, i, cnt, p: integer;
  s, temp: string;
  r: real;
  stop: boolean;
  name: array[1..cond]of string;                  // holds file names
  val: array[1..cond, 1..n, 1..k]of real;         // holds raw input data
  maxi: array[1..cond, 1..n]of integer;
  maxj, total: array[1..cond]of integer;
  curve: array[1..cond, 0..entries]of real;       // holds calculated kernel density curve values
  paired: array[1..cond, 1..n, 1..pairs]of real;  // holds all event pairs
  sumpair: array[1..cond, 1..n]of real;
  countpair: array[1..cond, 1..n]of integer;
  pairtemp: array[1..cond, 0..n * pairs]of real; 

begin
  // format arrays
  for y := 1 to cond do
    for j := 1 to n do
      for i := 1 to k do
        val[y, j, i] := -1;
  for y := 1 to cond do
    for j := 1 to n do
      maxi[y, j] := 0;
  for y := 1 to cond do
    maxj[y] := 0;
  for y := 1 to cond do
    total[y] := 0;
  for y := 1 to cond do
    for j := 1 to n do
      for i := 1 to pairs do
        paired[y, j, i] := -1;
  for y := 1 to cond do
    for j := 1 to n do
    begin
      sumpair[y, j] := 0;
      countpair[y, j] := 0;
    end;
  
  // set input file names
  name[1] := '52con.csv';
  name[2] := '52fluox.csv';
  name[3] := '72con.csv';
  name[4] := '72fluox.csv';
  name[5] := '96con.csv';
  name[6] := '96fluox.csv';
  
  // import data
  for y := 1 to cond do
  begin
    assign(fin, name[y]);
    reset(fin);
    j := 0;
    while(not eof(fin)) do
    begin
      readln(fin, s);
      inc(j);
      stop := FALSE;
      i := 0;
      while((length(s) > 0) and (not stop)) do
      begin
        if(s[1] = ',') then
          stop := TRUE
        else
        begin
          if(pos(',', s) > 0) then
            temp := copy(s, 1, pos(',', s) - 1)
          else 
          begin
            temp := s;
            stop := TRUE
          end;
          inc(i);
          val[y, j, i] := StrToFloat(temp);
        end;
        delete(s, 1, pos(',', s));
      end;
      maxi[y, j] := i;
      total[y] := total[y] + i;
    end;
    maxj[y] := j;
    close(fin);
  end;
  
  // generate kernel density curves
  for y := 1 to cond do
  begin
    r := rmin;
    while(r <= rmax) do
    begin
      cnt := 0;
      for j := 1 to maxj[y] do
        for i := 1 to maxi[y, j] do
          if((val[y, j, i] >= (r - window)) and (val[y, j, i] <= (r + window))) then
            inc(cnt);
      curve[y, round(r / fine)] := (cnt * fine) / (total[y] * 2 * window);
      r := r + fine;
    end;
  end;
  
  // set kernel density otput file names
  name[1] := '52out.csv';
  name[2] := '72out.csv';
  name[3] := '96out.csv';
  
  // write kernel density outputs
  for y := 1 to round(cond / 2) do
  begin
    assign(fout, name[y]);
    rewrite(fout);
    writeln(fout, '"TZ%","Con","Fluox"');
    r := rmin;
    while(r <= rmax) do
    begin
      writeln(fout, r, ',', curve[2 * y - 1, round(r / fine)], ',', curve[2 * y, round(r / fine)]);
      r := r + fine;
    end;
    close(fout);
  end;
  
  // compute pairs
  for y := 1 to cond do
    for j := 1 to maxj[y] do
      if maxi[y, j] > 1 then
      begin
        cnt := 0;
        for i := 2 to maxi[y, j] do
          for p := 1 to (i - 1) do
          begin
            inc(cnt);
            paired[y, j, cnt] := abs(val[y, j, i] - val[y, j, p]);
            sumpair[y, j] := sumpair[y, j] + paired[y, j, cnt];
          end;
        countpair[y, j] := cnt;
      end;
  
  // format pair distances for outputting
  for y := 1 to cond do
  begin
    cnt := 0;
    for j := 1 to maxj[y] do
      if(countpair[y, j] > 0) then
        for i := 1 to countpair[y, j] do
        begin
          inc(cnt);
          pairtemp[y, cnt] := paired[y, j, i];
        end;
    pairtemp[y, 0] := cnt;
  end;
  
  // write pair distance outputs
  cnt := 0;
  for y := 1 to cond do
    if(pairtemp[y, 0] > cnt) then
      cnt := round(pairtemp[y, 0]);
  
  assign(fout, 'pairdist.csv');
  rewrite(fout);
  writeln(fout, '"52Con","52Fluox","72Con","72Fluox","96Con","96Fluox"');
  
  for i := 1 to cnt do
  begin
    temp := '';
    for y := 1 to cond do
    begin
      if(i <= pairtemp[y, 0]) then
        temp := temp + FloatToStr(pairtemp[y, i]);
      temp := temp + ',';
    end;
    writeln(fout, temp);
  end;
  close(fout);
  
  // write average pair distance output
  cnt := 0;
  for y := 1 to cond do
    if(maxj[y] > cnt) then
      cnt := maxj[y];
  
  assign(fout, 'pairaverage.csv');
  rewrite(fout);
  writeln(fout, '"52Con","52Fluox","72Con","72Fluox","96Con","96Fluox"');
  
  for j := 1 to cnt do
  begin
    temp := '';
    for y := 1 to cond do
    begin
      if((j <= maxj[y]) and (countpair[y, j] > 0)) then
        temp := temp + FloatToStr(sumpair[y, j] / countpair[y, j]);
      temp := temp + ',';
    end;
    writeln(fout, temp);
  end;
  close(fout);
  
end.