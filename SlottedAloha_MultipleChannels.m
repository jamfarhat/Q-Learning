function [ntput, ndist, ntotal] = SlottedAloha_MultipleChannels(Devices, Relays, Channels_Relays, runs, frames, Slots, SNR, N, r)
% Loop Frames
ThroughputRuns=[];
TotalTraffic=zeros(1, runs);
TotalTrafficDistinct=zeros(1, runs);

for(l=1:frames)
    SuccessTransmission=zeros(Devices, runs);
    ThroughputFrame=zeros(1, runs);
    SlotChoosen=randi(Slots, Devices, runs);
    ChannelChoosen=randi(Channels_Relays, Devices, runs);

    for(k=1:Slots)
        for(s=1:runs)
            TransmittingDevices=find(SlotChoosen(:, s)==k);
            TransmittingChannel=ChannelChoosen(TransmittingDevices, s);

            if(length(TransmittingDevices)>=1)
                SNR_Device=SNR(TransmittingDevices', :, s);
                uniqueChannels=unique(TransmittingChannel);

                for(c=1:length(uniqueChannels))
                    SNR_Device_Channel=SNR_Device(TransmittingChannel'==uniqueChannels(c), :);
                    TransmittingDevices_channel=TransmittingDevices(TransmittingChannel'==uniqueChannels(c), :);
                    [SNR_Device_ord indexes]=sort(SNR_Device_Channel, 1, 'descend');
                    TransmittingDevices_ord=TransmittingDevices_channel(indexes);

                    uniqueRelays=1:1:Relays;

                    for(rr=1:length(uniqueRelays))
                        SIC_boolean=0;
                        for(jj=1:size(SNR_Device_Channel,1))
                            Interference=sum(SNR_Device_ord((jj+1):end, uniqueRelays(rr)));
                            SINR=(SNR_Device_ord(jj, uniqueRelays(rr)))./(Interference + N);

                            if(log2(1+SINR)>=r && SIC_boolean==0)
                                ThroughputFrame(s)=ThroughputFrame(s)+1;
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
    
    ThroughputRuns=[ThroughputRuns; mean(sum(SuccessTransmission>0))/Slots];
    TotalTrafficDistinct=[TotalTrafficDistinct+sum(SuccessTransmission>0)];
    TotalTraffic=[TotalTraffic+sum(SuccessTransmission)];
end

ndist=mean(TotalTrafficDistinct)/frames;
ntotal=mean(TotalTraffic)/frames;
ntput=(r/Channels_Relays)*mean(ThroughputRuns)*ndist/ntotal;