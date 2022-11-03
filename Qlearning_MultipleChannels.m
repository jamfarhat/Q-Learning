function [ntput, ndist, ntotal] = Qlearning_MultipleChannels(Devices, Relays, Channels_Relays, runs, frames, Slots, SNR, N, r, QTable, alpha, gamma)
% Loop Frames
ThroughputRuns=[];
TotalTraffic=zeros(1, runs);
TotalTrafficDistinct=zeros(1, runs);

for(l=1:frames)
    ThroughputFrame=zeros(1, runs);
    SuccessTransmission=zeros(Devices, runs);
    %Qtable
    DeviceTransmitting=zeros(1, Devices);
    %Q-learning Search
    for(dd=1:Devices)
        for(rr=1:runs)
            QTableDevice=QTable(:,:,rr,dd);
            maximum = max(max(QTableDevice));

            [y,x]=find(QTableDevice==maximum);

            if(length(x)>1)
                randomChoice=randi(length(x));

                ChannelChoosen(dd, rr)=y(randomChoice);
                SlotChoosen(dd, rr)=x(randomChoice);
                MaxQ(dd, rr)=maximum;
            else
                ChannelChoosen(dd, rr)=y;
                SlotChoosen(dd, rr)=x;
                MaxQ(dd, rr)=maximum;
            end
        end
    end

    Reward=-ones(Devices, runs);

    for(k=1:Slots)
        for(s=1:runs)
            TransmittingDevices=find(SlotChoosen(:, s)==k);
            TransmittingChannel=ChannelChoosen(TransmittingDevices, s);

            if(length(TransmittingDevices)>=1)
                SNR_Device=SNR(TransmittingDevices', :, s);
                uniqueChannels=unique(TransmittingChannel);

                for(c=1:length(uniqueChannels))
                    SNR_Device_Channel=SNR_Device(TransmittingChannel'==uniqueChannels(c), :);
                    [SNR_Device_ord indexes]=sort(SNR_Device_Channel, 1, 'descend');
                    TransmittingDevices_Channel=TransmittingDevices(TransmittingChannel'==uniqueChannels(c));

                    TransmittingDevices_ord=TransmittingDevices_Channel(indexes);

                    uniqueRelays=1:1:Relays;

                    for(rr=1:length(uniqueRelays))
                        SIC_boolean=0;
                        for(jj=1:size(SNR_Device_Channel,1))
                            Interference=sum(SNR_Device_ord((jj+1):end, uniqueRelays(rr)));
                            SINR=(SNR_Device_ord(jj, uniqueRelays(rr)))./(Interference + N);

                            if(log2(1+SINR)>=r && SIC_boolean==0)
                                ThroughputFrame(s)=ThroughputFrame(s)+1;
                                Reward(TransmittingDevices_ord(jj, uniqueRelays(rr)), s)=1;
                                SuccessTransmission(TransmittingDevices_ord(jj, uniqueRelays(rr)), s)=SuccessTransmission(TransmittingDevices_ord(jj, uniqueRelays(rr)), s)+1;
                            else
                                SIC_boolean=1;
                            end
                        end
                    end
                end
            end
        end
    end

    for(dd=1:Devices)
        for(rr=1:runs)
            QTableDevice=QTable(ChannelChoosen(dd, rr), SlotChoosen(dd, rr), rr, dd);
            ValorQTable=(1-alpha)*QTableDevice+alpha*(Reward(dd, rr)+gamma*QTableDevice);
            QTable(ChannelChoosen(dd, rr), SlotChoosen(dd, rr), rr, dd)=ValorQTable;
        end
    end

    ThroughputRuns=[ThroughputRuns; mean(sum(SuccessTransmission>0))/Slots];
    TotalTrafficDistinct=[TotalTrafficDistinct+sum(SuccessTransmission>0)];
    TotalTraffic=[TotalTraffic+sum(SuccessTransmission)];
end

ndist=mean(TotalTrafficDistinct)/frames;
ntotal=mean(TotalTraffic)/frames;
ntput=(r/Channels_Relays)*mean(ThroughputRuns)*ndist/ntotal;