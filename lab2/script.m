clc;


function res=f1(x, y)
  res = 1/(x + y^2);
endfunction

function res=f2(x, y)
  res = (y^4 + y.*sin(x))./cos(x);
endfunction

function res=f3(x, y)
  res = (x^2 + y^3)/(x*y^2);
endfunction

function res=rungeKutta(a, b, h, x0, y0, funcname)
  res(1) = y0;
  j = 2;
  for i=a:h:b - h
    k1 = feval(funcname, x0, y0);
    k2 = feval(funcname, x0 + h/2, y0 + (h/2)*k1);
    k3 = feval(funcname, x0 + h/2, y0 + (h/2)*k2);
    k4 = feval(funcname, x0 + h, y0 + h*k3);
    dy = (h/6) * (k1 + 2*k2 + 2*k3 + k4);
    y0 = y0 + dy;
    x0 = x0 + h;
    res(j) = y0;
    j += 1;
  endfor
endfunction


function res=kuttaMerson(a, b, step, x0, y0, eps, funcname)
  res(1) = y0;
  j = 2;
  for i=a:step:b - step
    h = step;
    while (1)
      k1 = feval(funcname, x0, y0);
      k2 = feval(funcname, x0 + h/3,y0 + (3/2)*k1);
      k3 = feval(funcname, x0 + h/3, y0 + (h/6)*k1 + (h/6)*k2);
      k4 = feval(funcname, x0 + h/2, y0 + (h/8)*k1 + ((3*h)/2)*k2);
      k5 = feval(funcname, x0 + h, y0 + (h/2)*k1 - ((3*h)/2)*k3 + 2*h*k4);
      nx = y0 + (h/2)*(k1 - 3*k3 + 4*k4);
      nx2 = y0 + (h/6)*(k1 + 4*k3 + k5);
      R = 0.2*abs(nx - nx2);
      if (R > eps)
        h = h/2;
      else
        res(j) = nx2;
        j += 1
        x0 += h;
        y0 = nx2;
        break;
      endif
    endwhile
  endfor
endfunction



y0 = 1;
span = [0 1];
[od1x, od1y] = ode23('f1', span, y0);
[odlx2, odly2] = ode45('f1', span, y0);

x = 0:0.1:1;
s1 = rungeKutta(0, 1, 0.1, 0, 1, 'f1');
s2 = kuttaMerson(0, 1, 0.1, 0, 1, 0.001, 'f1')

plot(x, s1, x, s2, od1x, od1y, odlx2, odly2);



























