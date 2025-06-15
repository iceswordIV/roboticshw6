function a = hatm_num(p)

a1 = p(1,1);
a2 = p(2,1);
a3 = p(3,1);

a = [
    0, -a3, a2;
    a3, 0, -a1;
    -a2, a1, 0
    ];

end
