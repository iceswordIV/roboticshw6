function Adg_ab = adjmginv(g_ab)
R_ab = [
    g_ab(1,1), g_ab(1,2), g_ab(1,3);
    g_ab(2,1), g_ab(2,2), g_ab(2,3);
    g_ab(3,1), g_ab(3,2), g_ab(3,3)
    ];

p_ab = [g_ab(1,4); g_ab(2,4); g_ab(3,4)];

a1 = p_ab(1,1);
a2 = p_ab(2,1);
a3 = p_ab(3,1);

phat = [
    0, -a3, a2;
    a3, 0, -a1;
    -a2, a1, 0
    ];



Adg_ab = [
    R_ab', -R_ab' * phat;
    zeros(3,3), R_ab'
    ];

end