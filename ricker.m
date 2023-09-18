


function y = ricker(f,t,to);

y = (1-(2.0*((pi)^2)*(f^2)*(t-to)^2))*(exp(-((pi)^2)*(f^2)*(t-to)^2));

end
