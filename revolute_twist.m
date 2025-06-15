function twist = revolute_twist(w,q)

twist = [cross(-w,q);w];
end
