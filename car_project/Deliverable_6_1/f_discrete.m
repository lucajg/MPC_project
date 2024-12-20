function xp = f_discrete(x,u,Ts,f)
    xp = RK4(x,u,Ts,f);
end