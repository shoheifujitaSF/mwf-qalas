function magnetization = pulse_ode45(magnetization, Bx,By,Gz, R2, R1, gamma, r, b0, b1, rho,time_grid)

for rr = 1:length(r)
    %-Defining ODE function
    A  = @(t) [-R2, gamma*Gz(t)*r(rr) + b0(rr), -gamma*By(t)*b1(rr), 0;...
        -(gamma*Gz(t)*r(rr) + b0(rr)), -R2, gamma*Bx(t)*b1(rr), 0;...
        gamma*By(t)*b1(rr), -gamma*Bx(t)*b1(rr), -R1, rho * R1;...
        0, 0, 0, 0];
    
    ode_func = @(t,y) A(t)*y;
    
    %-ODE 45
    [~,y] = ode45(ode_func,time_grid,magnetization(:,rr));

    magnetization(:,rr) = squeeze(y(end,:));
end

end