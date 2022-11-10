function r = rmse( x_es , x ) 
    r = 0 ;
    for i = 1 : length(x) 
       r = r + (x_es(i) - x(i) ).^2 ;
    end
    r = sqrt(r/length(x));
end