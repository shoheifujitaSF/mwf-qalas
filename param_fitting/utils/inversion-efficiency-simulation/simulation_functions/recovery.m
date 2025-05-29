function magnetization = recovery(magnetization,R2,R1,rho,time)

    magnetization(1:2,:) = magnetization(1:2,:) * exp(-time * R2);
    magnetization(3,:)   = magnetization(3,:) * exp(-time * R1) + rho * (1 - exp(-time * R1));
end