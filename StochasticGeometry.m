function c = StochasticGeometry(Devices, Relays, Radius, runs)
    for jj=1:Relays
        %----------------------------------------------------
        Raio = Radius*sqrt(rand(1,runs));
        Theta = 2*pi*rand(1,runs);
        %----------------------------------------------------
        positionsX = Raio.*cos(Theta);
        positionsY = Raio.*sin(Theta);
        %----------------------------------------------------
        relay_position_x(jj, :)=positionsX;
        relay_position_y(jj, :)=positionsY;
    end
 
    for ii=1:Devices
        %----------------------------------------------------
        Raio = Radius*sqrt(rand(1,runs));
        Theta = 2*pi*rand(1,runs);
        %----------------------------------------------------
        positionsX = Raio.*cos(Theta);
        positionsY = Raio.*sin(Theta);
        %----------------------------------------------------
        for(jj=1:Relays)
            c(ii, jj, :)=sqrt((positionsX-relay_position_x(jj, :)).^2 + (positionsY-relay_position_y(jj, :)).^2);
        end
    end
end