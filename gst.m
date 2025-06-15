function gst = gst(q)



I = [
    1,0,0;
    0,1,0;
    0,0,1
    ];

gst = [
    I q
    zeros(1,3) 1
    ];

end