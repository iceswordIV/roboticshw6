function gsl_0 = gsl_0(center_point)



I = [
    1,0,0;
    0,1,0;
    0,0,1
    ];

gsl_0 = [
    I center_point
    zeros(1,3) 1
    ];

end