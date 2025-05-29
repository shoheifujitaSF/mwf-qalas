function magnetization = t2prep(magnetization, event_blocks,rho,R1,R2,b0,b1,r,gamma)

for ii = 1:length(event_blocks)
    if(event_blocks{ii}.pulseflag == 1)
        %-simulate a pulse with ode 45
        Bx = event_blocks{ii}.Bx;
        By = event_blocks{ii}.By;
        Gz = event_blocks{ii}.Gz;
        time_grid = event_blocks{ii}.time_grid;

        magnetization = pulse_ode45(magnetization, Bx,By,Gz, R2, R1, gamma, r, b0, b1, rho, time_grid);
    else
        time = event_blocks{ii}.time_grid(end) - event_blocks{ii}.time_grid(1);

        %-simulate recovery
        magnetization = recovery(magnetization,R2,R1,rho,time);
    end
end
end