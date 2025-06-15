function etheta = exponential_epison_hat_theta(w,v,theta)

I = [
    1,0,0;
    0,1,0;
    0,0,1
    ];


etheta = [
    exponential(w,theta), (I - exponential(w,theta)) * cross(w,v) ;
    zeros(1,3), 1
    ];

end
