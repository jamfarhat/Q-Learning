function c = InitializeQTable(Devices, Relays, Slots, Runs, Initialization)
    if(Initialization==true)
        c=zeros(Relays, Slots, Runs, Devices);
    else
        c=-1+(2)*rand(Relays, Slots, Runs, Devices);
    end
end