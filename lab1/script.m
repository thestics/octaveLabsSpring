# 123
clc;

function res = f(x, y)
  res = (y*x^2 - x^2)/(y^2 + x*y^2);
endfunction

 
function res = solveDiffEq(a, b, h, x0, y0) # euler
  j = 1;
  for i=a:h:b
    x1 = x0 + h; 
    y1 = y0 + h*f(x0, y0);
    res(j) = y1;
    j += 1;
    x0 = x1;
    y0 = y1;
  endfor
endfunction


function res = solveDiffEq2(a, b, h, x0, y0)  # hoit
  j = 1;
  for i=a:h:b
    x1 = x0 + h; 
    y1 = y0 + h*(f(x0, y0) + f(x0 + h, y0 + f(x0, y0)));
    res(j) = y1;
    j += 1;
    x0 = x1;
    y0 = y1;
  endfor
endfunction


function main()
  s = solveDiffEq(1, 4, 0.1, 1, 2);
  s1 = solveDiffEq2(1, 4, 0.1, 1, 2);
  s2 = solveDiffEq(1, 4, 0.05, 1, 2);
  x = 1:0.1:4;
  x2 = 1:0.05:4;
  plot(x, s, x, s1, x2, s2);
endfunction
