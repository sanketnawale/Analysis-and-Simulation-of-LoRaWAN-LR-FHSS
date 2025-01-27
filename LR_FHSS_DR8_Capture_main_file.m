% Paper title: Analysis and Simulation of LoRaWAN LR-FHSS for Direct-to-Satellite Scenario
% IEEE XPlore: https://ieeexplore.ieee.org/document/9653679
% Authors: Muhammad Asad Ullah, Konstantin Mikhaylov, Hirley Alves

% Cite this: M. A. Ullah, K. Mikhaylov and H. Alves, "Analysis and Simulation of LoRaWAN LR-FHSS for Direct-to-Satellite Scenario," in IEEE Wireless Communications Letters, doi: 10.1109/LWC.2021.3135984.

clear all 
close all

%% Packet and transmissions parameters
tic
Payload = 10;           % Message payload
Header_N_DR8 = 3;       % Header replicas
Code_Rate = 1/3;
Header_duration = 0.233; %233 ms long headers
F_duration = 0.05;       %50 ms payload data fragments
Header_ToA_DR8 = Header_N_DR8*Header_duration;
%Nodes = [50 60 70 80 90 100].*1e3;
Nodes = 50e3;
Simulation_T = 3600; % 1 hour duration
pkct_p_h = 4;      % Packets per hour per end-device
OBW_channels=280;  % No. of OBW channels
MonteCarlo = 1e3;    % No. of Iterations
M = 2;             % M = 2 for DR8, M=4 for DR9

%% Gains and Pt are converted into linear form

Pt = 10^(14/10)/1000;      % Transmit Power of LoRa 14 dBm
Freq_Band = 868e6;         % 868 MHz (frequency band Europe)
SpeedLight  = 3e8;         % Speed of light
wavelength = SpeedLight/Freq_Band;

Gr=(10.^((22.6)/10));      %22.6 dBi: Gateway at Satellite
Gt=(10.^((2.15)/10));      %2.15 dBi: End-device
eta = 2;

%% parameters for satellites
E = 10:1:90;               %Elevation Angles
R = 6378e3;                % Radius of earth
H = 780e3;                 %Orbital height  

%% Distance from user to satellite as function of elevation angle
[Distance, Elevation_Angles, Ground_distance,FootPrint_R]=Satellite_Geometry(H,E);
%distance from user to satt as a function of elevation angle 
%E: ele sangle 10:1:90
%dist: dist in staraight bline from tghe user to the sat at each elevation angle
%elevation : ele angle are analyzed
%gnd distance1: how far hori ythe user signal trave; on eath surface
%foothprintyR: the radius of the sat coverage foothprint on earth
X = [1 8 11 14 17 19 22 25]; %To simulate fewer points : specific elevation angles reducingthe dataset
Distance=Distance(X); 
% filters the bdistance array to include  only the distance 
% corresponding to these selected indices.

E_angles = [10 20 30 40 50 60 70 80 90];
%an anglgls taht how sat is high in sky
K_factor = [1.24 3.07 3.24 3.6 3.89 5.63 9.77 17.06 25.11];
%k: describe teh fadinga due to atmosphere and geometry path
k = sort(interp1(E_angles,K_factor,Elevation_Angles),'descend'); % interpolation of fineer granulatitu
%interp1: this interpolation esimate the k factor for every Eangle
%ip : known E angle and kfactor
%op interpolated kfactor for all anglks in elevation angles from sat geo calculation
%sort:sort the interpolated kfactors in decending order

%% Time on air function calculte the toa for tx packt in lorawan
[ToA_DR8,ToA_DR8_WH] = ToA_Packets_DR8(Payload,Header_ToA_DR8,M); 
%ipo payload in bytesa
%header are precalculated , M is spreafding factor
% ToA_DR8 -> including headers duration
% ToA_DR8_WH -> without headers duration
Transceiver_wait = 6.472/1000; %Tw      dealy by transreviver known asd tw
ToA_DR8(1)=ToA_DR8(1) + Transceiver_wait;  % Total on-air time

%% Number of fragments each fragmenr small data pkt has fixed duration 50nms
fragment_duration = 50/1000; %50 ms

fragment_50_ms = floor(ToA_DR8_WH(1)/fragment_duration);  %No. of payload data fragments
% The last fragment may not be equal to 50ms. We need to calculate it.
% determine how many full fragmnet 50ms fit intro teh total payload time  on air (TOA dr8WH
Last_fragment_duration = ((ToA_DR8_WH(1)/fragment_duration) - fragment_50_ms)*fragment_duration;
%calculte teh leftover duration for thr last fragment which dont fit in 50ms
fragment_PHY_length  = fragment_50_ms + 1;
%calcultr th total no of frangments + accounts for the last fraegmnets
fragment_length = Header_N_DR8 + length(Transceiver_wait) + fragment_PHY_length;
%total frangment length
%pack_tx_segments = Header_N_DR8 + ceil((0.102*ceil((Payload+2)/4))/F_duration) +length(Transceiver_wait);

%% Simulator
transmission_sim=zeros(MonteCarlo,fragment_length);%initilization
%montecarlo: no of iteration : fragment length : total no of fragments we calulate eralier

%llop over distances
for c=1:1:length(Distance) %C index for each distance in the distance arry
N=Nodes;
csvwrite('loop_count.txt',c)
decoded = 0;% Number of successfully decoded packets
decoded_Capture = 0; % Successfully decoded packets with capture effect
discarded = 0;% Packets discarded due to collisions
H_success = 0; % Successful header transmissions
F_success = 0;% Successful fragment transmissions
    for m=1:1:MonteCarlo 
    %inner loop for monte carlo
     clear TX_stamp
     clear Tsort
     clear TimeStamp
     clear pack_tx_segments

     First_transmit = rand(1,1)/1000;          % First transmission time in sec
     %Mean Inter-arrival Time       
     mu = (1/(N*pkct_p_h)).*Simulation_T;      % inter arrival time
     Inter_arrivals = exprnd(mu,1,N*pkct_p_h); % inter arrival of the traffic in a hour
     %Generate random inter-arrival times using an exponential distribution
     Times = [First_transmit Inter_arrivals];%Convert inter-arrival times into cumulative timestamps for when each packet is transmitted.
     
     %Next step: convert inter-arrivals into time stamp 
     TimeStamp = cumsum(Times);                % Time stamp of the traffic in the network
     %cumsum(Times): Calculates the cumulative sum to create the timestamps
     pack_tx_segments=zeros(length(TimeStamp),fragment_length);%Packet Transmission Segments
     %Create a matrix to store the timestamp for each fragment of every packet.
     %%Rows represent packets.
     %Columns represent fragments (headers + payload fragments).
     
     %% Time stamp of the hops (segments)
         %Iterating Over Packets   
        for pack=1:1:length(TimeStamp)   %For each packet, calculate the timestamps for its fragments.
            for frag = 1:1:fragment_length %Calculate the timestamp for each fragment of a packet based
                if frag == 1%The first fragment's timestamp is simply the packet's start tim
                pack_tx_segments(pack,frag) = TimeStamp (pack);%The first fragment's timestamp is simply the packet's start tim
                elseif frag > 1 && frag <=(Header_N_DR8+1)%Add Header_duration (fixed at 0.233 seconds) to the previous fragment’s timestamp
                pack_tx_segments(pack,frag) = pack_tx_segments(pack,frag-1) + Header_duration;
                elseif frag == (Header_N_DR8+2)%After the headers, there is a small wait time (6.472 ms) before payload fragments are transmitted.
                pack_tx_segments(pack,frag) = pack_tx_segments(pack,frag-1) + Transceiver_wait;
                else %The payload fragments are transmitted one by one, with a fixed duration of fragment_duration (50 ms) for each fragment
                pack_tx_segments(pack,frag) = pack_tx_segments(pack,frag-1) + fragment_duration;
                end
            end
        end

    %% Vulnerable time
    Transmit = randi(length(TimeStamp));  % Randomly selects one of the devices that has transmitted.
    % Ts: The time when this selected device starts transmitting.
    Ts= TimeStamp(Transmit);
    
    % vulnerable time
    Tstart = Ts - ToA_DR8(1);%Time when interference can start. This is the transmission start time minus the packet's Time on Air (ToA_DR8(1)).
    Tend = Ts + ToA_DR8(1);% Time when interference can end. This is the transmission end time plus the Time on Air.

%% Find the number active devices when the desired device was transmitting

    index_start = find(pack_tx_segments(:,1)>=Tstart & pack_tx_segments(:,1)<=Tend); % Select all the transmission in 2T interval, where T is on-air time
%pack_tx_segments(:,1): The first column of pack_tx_segments contains the start times of all transmissions.
    %find(...): Selects transmissions that started between Tstart and Tend
    simultaneous = (unique([index_start'])); %simultaneous: A list of all devices transmitting during the "vulnerable time.
    target_index = find(simultaneous == Transmit);%target_index: The index of the selected device in the simultaneous list.

    simultaneous(target_index) = [];%Remove the selected device to avoid comparing it with itself.
 
%% Frequency-time scheduling of target transmission: which fragment is using the specific channel for a specific time?
% Initialize the schedule
    target_pattern=zeros(1,size(pack_tx_segments,2));%An array representing which frequency channel is used for each fragment.
    target_pattern(1) = randi(OBW_channels,1,1);            %First hop of the desired signal.. Assign a random channel for the first fragment.

        for assign=2:1:size(pack_tx_segments,2)%Assign a frequency channel for each subsequent fragment.
            if(assign==Header_N_DR8+1)     %          
                target_pattern(assign) = 0;                % Do not assign a channel during T_wait :::During the Transceiver Wait (after the headers), no frequency channel is assigned.
            else  
                target_pattern(assign) = randi(OBW_channels,1,1);%Assign random channels for other fragm
                dif_track=abs(target_pattern(assign)-target_pattern(assign-1));
                while dif_track < 8                       % 8 x 488 = 3.9 kHz spacing ;;If the difference (dif_track) is less than 8, the channel spacing is too small, and interference could occur.
                    target_pattern(assign) = randi(OBW_channels,1,1);
                    dif_track=abs(target_pattern(assign)-target_pattern(assign-1));
                end
            end 
        end
    
    %% Collision analysis
    target_collided = zeros(1,size(pack_tx_segments,2));           %Counts how many times each fragment (header or payload) collided.
    target_discarded = zeros(1,size(pack_tx_segments,2));          %Tracks fragments that can't be used even with recovery techniques (e.g., capture effect).

    frag_active = zeros(1,length(pack_tx_segments(Transmit,:)));%Tracks which fragments are active (being transmitted) at each moment.

    
    clear seg_simultaneous %Extracts the transmission segments of all devices transmitting during the vulnerable time.
    seg_trans_simultaneous = pack_tx_segments(simultaneous,:);                % Select all the fragments during the interval 2T
    seg_trans_simultaneous(:,(Header_N_DR8+1))=[];                            %Twait, that's not a transmission : %Removes the Transceiver Wait (T_wait) fragment since it doesn’t involve any actual transmissi
    seg_simultaneous = seg_trans_simultaneous;
    %% Following Equation (4), (5), (6) to find A_{H}, A_{F} and A_{L}
        for seg=1:1:length(pack_tx_segments(Transmit,:))%Go through each fragment of the selected device (both headers and payloads).
            
            if(seg~=Header_N_DR8+1) % Ignore the Transceiver Wait (T_wait) fragment because it’s not part of an active transmission.
                clear active_pattern
                clear iscollision

                if(seg~=length(pack_tx_segments(Transmit,:)) && seg<=Header_N_DR8)         % collisions during header transmission
                    
                transmission_sim_header(m,seg) = length(find(seg_simultaneous(:,1:Header_N_DR8)>=(pack_tx_segments(Transmit,seg)-Header_duration) & seg_simultaneous(:,1:Header_N_DR8)<=(pack_tx_segments(Transmit,seg)+Header_duration))); 
                %Checks how many other devices' headers overlap in time with the current header.Finds overlaps where other headers started within the current header's time range.
                
                transmission_sim_frag(m,seg) = length(find(seg_simultaneous(:,(Header_N_DR8+1):end-1)>=(pack_tx_segments(Transmit,seg)-F_duration) & seg_simultaneous(:,(Header_N_DR8+1):end-1)<=(pack_tx_segments(Transmit,seg)+Header_duration)));
                %Checks if any payload fragments from other devices overlap with the current header.
                transmission_sim_last(m,seg) = length(find(seg_simultaneous(:,end)>=(pack_tx_segments(Transmit,seg)-Last_fragment_duration) & seg_simultaneous(:,end)<=(pack_tx_segments(Transmit,seg)+Header_duration)));
                %Checks if the last fragment of other devices overlaps with the current header.
                transmission_sim(m,seg) = transmission_sim_header(m,seg) + transmission_sim_frag(m,seg) + transmission_sim_last(m,seg);
                %
                elseif(seg~=length(pack_tx_segments(Transmit,:)) && seg>Header_N_DR8)    % collisions during payload data fragments transmission
                 %Similar logic as for headers, but now focused on payload fragments:
   
                transmission_sim_header(m,seg) = length(find(seg_simultaneous(:,1:Header_N_DR8)>=(pack_tx_segments(Transmit,seg)-Header_duration) & seg_simultaneous(:,1:Header_N_DR8)<=(pack_tx_segments(Transmit,seg)+F_duration)));
                transmission_sim_frag(m,seg) = length(find(seg_simultaneous(:,(Header_N_DR8+1):end-1)>=(pack_tx_segments(Transmit,seg)-F_duration) & seg_simultaneous(:,(Header_N_DR8+1):end-1)<=(pack_tx_segments(Transmit,seg)+F_duration)));
                transmission_sim_last(m,seg) = length(find(seg_simultaneous(:,end)>=(pack_tx_segments(Transmit,seg)-Last_fragment_duration) & seg_simultaneous(:,end)<=(pack_tx_segments(Transmit,seg)+F_duration)));
                
                transmission_sim(m,seg) = transmission_sim_header(m,seg) + transmission_sim_frag(m,seg) + transmission_sim_last(m,seg);
                
                else                                                                     % collisions during the last fragment transmission
                transmission_sim_header(m,seg) = length(find(seg_simultaneous(:,1:Header_N_DR8)>=(pack_tx_segments(Transmit,seg)-Header_duration) & seg_simultaneous(:,1:Header_N_DR8)<=(pack_tx_segments(Transmit,end))+Last_fragment_duration));
                transmission_sim_frag(m,seg) = length(find(seg_simultaneous(:,(Header_N_DR8+1):end-1)>=(pack_tx_segments(Transmit,seg)-F_duration) & seg_simultaneous(:,(Header_N_DR8+1):end-1)<=(pack_tx_segments(Transmit,end))+Last_fragment_duration));
                transmission_sim_last(m,seg) = length(find(seg_simultaneous(:,end)>=pack_tx_segments(Transmit,seg) & seg_simultaneous(:,end)<=(pack_tx_segments(Transmit,end))+Last_fragment_duration));
                
                transmission_sim(m,seg) = transmission_sim_header(m,seg) + transmission_sim_frag(m,seg) + transmission_sim_last(m,seg);
                                
                end   
                if (~isempty(transmission_sim(m,seg))):Tells how many other devices are transmitting at the same time as the current fragment.
                %If it’s not empty, there are devices interfering.

                frag_active(seg) =  (transmission_sim(m,seg));%
                % Select the channels for the simultenous fragments
                active_pattern = randi(OBW_channels,(transmission_sim(m,seg)),1)';
                % How many active devices are assigned to the same channel (collision)
                iscollision = find(active_pattern==target_pattern(seg));    % if non zero, that's collisionIf any interfering device uses the same channel as the selected device, a collision occurs.
                    
                    if(~isempty(iscollision))
                      target_collided (seg) = 1; %Marks the current fragment (seg) as "collided" because there’s interference on the same channe
                            
                      clear Coordinates                    
                      clear rho
                      clear Theta
                      %Generate random positions for interfering devices
                      Coordinates=zeros(length(iscollision),2);  
                      
                      rho = sqrt(rand(length(iscollision),1).*(max(Ground_distance)^2)); %Random distance of the device from the center of the coverage area.       
                      Theta = rand(length(iscollision),1)*2*pi; %     Random angle to place the device in a circular area                
                      Coordinates(:,1) = cos(Theta).*rho;% %Stores the (x, y) positions of the interfering devices.
                      Coordinates(:,2) = sin(Theta).*rho;
                      
                %% Generating location of interfering signals based on uniform distribution of NODES

                % (i)First location -> (ii) elevation angle -> (iii) Rician K factor
                %Calculate distances to the satellite
                      Location_Nodes_Int = sqrt(Coordinates(:,1).^2 + Coordinates(:,2).^2)'; %Distance of each interfering device from the center of the coverage area (on the ground)
                      dPropogation=zeros(1,length(iscollision)); %Total distance from the interfering device to the satellite, including the height of the satellite (H = 780 km)
                       E_dpro=zeros(1,length(iscollision)); % Intermediate step in calculating the angle using satellite geometry.
                       %Calculate elevation angles
                        % Distance from interfering nodes to satellite
                        for track=1:length(Location_Nodes_Int)
                        dPropogation(1,track) = sqrt(H^2 + Location_Nodes_Int(track).^2); 
                        end
                
                        for Fin_ang=1:length(Location_Nodes_Int)
                        E_dpro(Fin_ang) = (H*((H+2*R)) - dPropogation(Fin_ang).^2)./(2.*dPropogation(Fin_ang).*R);
                        end

                        E_AngPro=asind(E_dpro);%Converts the result to degrees using asind.

                %Determine the K-factor
                %Looks up the Rician K-factor for each interfering device based on its elevation angle (E_AngPro).
%K_factor: Represents how much interference can be mitigated due to line-of-sight conditions.
%interp1: Interpolates the K-factor for the specific elevation angle.
                        kC = interp1(E_angles,K_factor,E_AngPro);
                
                        %% Rician fading for interfering signals
               
                        muC = sqrt(kC./(2.*(kC+1)));    % Mean 
                        sigmaC = sqrt(1./(2.*(kC+1)));  % Variance 
                        hrC=(sigmaC.*randn(1,length(iscollision))) + muC;
                        hiC=1j.*(sigmaC.*randn(1,length(iscollision)) + muC);
                
                        h1C=(abs(hrC+hiC)).^2;
                
                        %% Interfering signal power

                        % Pr = pt*Gt*Gr*path loss * h  
    
                        %% total interferance = sum of all the interfering signals 
                        pr_h_g_I = sum(Pt.*h1C.*Gr.*Gt.*((wavelength./(4*pi.*dPropogation)).^eta));
                
                        %% Rician fading for desired signals
                
                        kD = k(c);
                
                        muD = sqrt(kD./(2*(kD+1)));     % Mean 
                        sigmaD = sqrt(1./(2*(kD+1)));   % Variance 
                
                        hrD=sigmaD*randn(1,1)+muD;
                        hiD=1j.*(sigmaD*randn(1,1)+muD);
                        h1D=(abs(hrD+hiD)).^2;
    
                        %% Received power of desired signal
               
                        pr_h_g_D = Pt.*h1D.*Gr.*Gt.*((wavelength./(4*pi.*Distance(c))).^eta);
                  
                       % https://lora-developers.semtech.com/library/tech-papers-and-guides/lora-and-lorawan/
                

                            if  pr_h_g_D  < (pr_h_g_I*4)
                            %discarded = discarded + 1; %% Destructive collision
                            target_discarded(seg) = 1;
                            end

                    end
                end
            end
        end


%% Decoding
% first three are header
% the fourth is t_wait
% Rest are fragments of 50 ms

Success_header = Header_N_DR8 - length(nonzeros(target_collided(1:Header_N_DR8)));       % No. of successfully received headers
Threshold = size(pack_tx_segments,2) - round(fragment_PHY_length *(1-Code_Rate))-Header_N_DR8 - length(Transceiver_wait);
Success_fragment = size(target_collided,2) - length(nonzeros(target_collided((Header_N_DR8+2):end)))-Header_N_DR8-1;

                      if(Success_fragment>=Threshold)
                           F_success=F_success+1;
                      end

        if (Success_header>=1)
            H_success=H_success+1;
            if(Success_fragment>=Threshold)
                decoded = 1 + decoded;                
            end
        end
        
Success_header_capture =Header_N_DR8 - length(nonzeros(target_discarded(1:Header_N_DR8)));

        if (Success_header_capture>=1)
        Success_fragment_capture = size(target_discarded,2) - length(nonzeros(target_discarded(Header_N_DR8+2:end)))-Header_N_DR8-1;
            if(Success_fragment_capture>=Threshold)
        decoded_Capture = 1 + decoded_Capture;
            end
    
        end        
        
     end
PS_DR8(c)=decoded;   %Simulated overall success probability
%PH_DR8(c,m)=H_success; %Simulated success probability of headers
%PF_DR8(c,m)=F_success; %Simulated success probability of data fragments

PS_DR8_Capture(c)=decoded_Capture; %Simulated success probability with capture effect

[PS_DR8_analytical,PH,PF] = DR8_analytical (N,pkct_p_h,Header_duration,F_duration,Last_fragment_duration,fragment_length,Header_N_DR8,Threshold,OBW_channels);
PA_S(c) = PS_DR8_analytical;  %Analytical overall success probability
%PA_H(c) = PH;                %Analytical success probability for headers
%PA_F(c) =PF;                  %Analytical success probability for data fragments
end


figure(2)
h(1)=plot(Distance/1e3,PA_S,'b-');
hold on
h(2)=plot(Distance/1e3,PS_DR8/MonteCarlo,'bo');
h(3)=plot(Distance/1e3,PS_DR8_Capture/MonteCarlo,'b--');
hold off
grid on
ylabel('Success probability', 'Interpreter', 'Latex');
xlabel('Distance from node to satellite (km)', 'Interpreter', 'Latex');
axis([Distance(1)/1e3 Distance(end)/1e3 0 1])
legend('DR8 analytical','DR8 Simulated', 'DR8 Capture effect')
toc
