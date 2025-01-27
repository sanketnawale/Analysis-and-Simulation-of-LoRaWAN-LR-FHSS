% Paper title: Analysis and Simulation of LoRaWAN LR-FHSS for Direct-to-Satellite Scenario
% IEEE XPlore: https://ieeexplore.ieee.org/document/9653679
% Authors: Muhammad Asad Ullah, Konstantin Mikhaylov, Hirley Alves

% Cite this: M. A. Ullah, K. Mikhaylov and H. Alves, "Analysis and Simulation of LoRaWAN LR-FHSS for Direct-to-Satellite Scenario," in IEEE Wireless Communications Letters, doi: 10.1109/LWC.2021.3135984.
%PS_DR8_analytical: The overall success probability (analytically) for Data Rate 8 (DR8).
%H_N_Pro_succ: Probability of successfully receiving all the header fragments
%P_F: Probability of successfully receiving enough payload fragments to decode the message.
function [PS_DR8_analytical,H_N_Pro_succ,P_F] = DR8_analytical (N,pkct_p_h,Header_duration,F_duration,Last_fragment_duration,fragment_length,Header_N_DR8,Threshold,OBW_channels)
%Inputs:
%N: Number of devices in the network.
%pkct_p_h: Packets sent per hour by each device.
%Header_duration: Time duration of each header fragment.
%F_duration: Time duration of each payload fragment.
%Last_fragment_duration: Time duration of the last fragment.
%fragment_length: Total number of fragments in the packet.
%%Header_N_DR8: Number of header fragments.
%Threshold: Minimum number of fragments needed for successful decoding.
%OBW_channels: Total number of frequency channels.

%% Probability of header success
% 1. Header collision with header
%Calculates how many header fragments from other devices overlap in time with a header fragment from the desired device.
H_active_2TH = ((N*pkct_p_h*(Header_N_DR8))/3600)*2*Header_duration;
% 2. Header collision with payload data fragments Collisions Between Header and Payload Fragments
H_active_TF_TH = ((N*pkct_p_h*(fragment_length-Header_N_DR8-1-1))/3600)*(Header_duration+F_duration);
%Calculates how many payload fragments from other devices collide with the header fragments of the desired device.
%Subtracts the number of headers, transceiver wait, and last fragment to count only payload fragments.
% 3. Header collision with last data fragmentCollisions Between Header and Last Fragment 
H_active_TL_TH = ((N*pkct_p_h)/3600)*(Header_duration+Last_fragment_duration);
%Calculates how many last fragments from other devices collide with the desired device’s header fragments
% Total colliding packet elements
% Total colliding packet elements
H_active = H_active_2TH + H_active_TF_TH + H_active_TL_TH;%Adds up all header collisions from the above steps.

%Header Success Probability
H_probability_succ=(1-(1/OBW_channels))^(H_active-1);%Probability that a single header fragment avoids collision, considering the number of channels.
H_N_Pro_succ = 1-(1-H_probability_succ).^Header_N_DR8;%Probability that at least one of the Header_N_DR8 header fragments avoids collision.

%% Fragments success probability, even if 1/3 (for DR8) of the fragments are lost

%% Payload data Fragment collisions
%1. Fragments collision with fragments
%Collisions Between Payload Fragments
F_active_12_TF = ((N*pkct_p_h*(fragment_length-Header_N_DR8-1-1))/3600)*2*F_duration;
%Calculates collisions between payload fragments from other devices and the desired device’s payload fragments.
%2. Fragments collision with headersCollisions Between Payload and Header Fragments
F_active_12_TH = ((N*pkct_p_h*(Header_N_DR8))/3600)*(Header_duration+F_duration);
%Calculates collisions between header fragments from other devices and the desired device’s payload fragments.
%3. Fragments collision with last data fragmentCollisions Between Payload and Last Fragment
F_active_12_TL = ((N*pkct_p_h)/3600)*(F_duration+Last_fragment_duration);
%Calculates collisions between last fragments from other devices and the desired device’s payload fragments.
% Total colliding packet elementsTotal Payload Collisions
F_active_12 = F_active_12_TF + F_active_12_TH + F_active_12_TL;%Adds up all payload fragment collisions and calculates the probability that a single payload fragment avoids collision.

%Last Fragment Success Probability


Frag_probability_succ=(1-(1/OBW_channels))^(F_active_12-1);           % excluding last fragment
%% Last fragment collisions
% 1. Last fragment collision with last fragments
F_active_last_TL = ((N*pkct_p_h)/3600)*2*(Last_fragment_duration);
% 2. Last fragment collision with headers
F_active_last_TH = ((N*pkct_p_h*(Header_N_DR8))/3600)*(Last_fragment_duration+Header_duration);
% 3. Last fragment collision with payload data fragments
F_active_last_TF = ((N*pkct_p_h*(fragment_length-Header_N_DR8-1-1))/3600)*(F_duration+Last_fragment_duration);
% Total colliding packet elements
F_active_last = F_active_last_TH +F_active_last_TF+F_active_last_TL;


%Calculates the probability that the last fragment avoids collision, considering collisions with:
%Other last fragments.
%Header fragments.
%Payload fragments.
Frag_probability_succ_last=(1-(1/OBW_channels))^(F_active_last-1);    %Last fragment

%Average Fragment Success

Total_fragment = fragment_length-Header_N_DR8-1;%Calculates the average success probability for all payload fragments, including the last one.
Frag_avg_suc = (sum(Frag_probability_succ*ones(1,(fragment_length-Header_N_DR8-2)))+Frag_probability_succ_last)/(Total_fragment);
%% Binomial probability
%Binomial Success Probability
Fragment_success = binopdf(1:1:Total_fragment,Total_fragment,Frag_avg_suc);
%binopdf(...): Models the probability of successfully receiving a certain number of fragments out of the total
%P_F: Probability of receiving enough fragments (greater than or equal to the Threshold)

%% Overall probability = Header success * fragment success
%Overall Success Probability
P_F = sum(Fragment_success(Threshold:1:Total_fragment));
PS_DR8_analytical = H_N_Pro_succ*sum(Fragment_success(Threshold:1:Total_fragment));
%Combines header and payload success probabilities to calculate the overall success probability for the message.
end
