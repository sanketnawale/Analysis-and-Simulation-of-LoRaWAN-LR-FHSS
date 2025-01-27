% Paper title: Analysis and Simulation of LoRaWAN LR-FHSS for Direct-to-Satellite Scenario
% IEEE XPlore: https://ieeexplore.ieee.org/document/9653679
% Authors: Muhammad Asad Ullah, Konstantin Mikhaylov, Hirley Alves

% Cite this: M. A. Ullah, K. Mikhaylov and H. Alves, "Analysis and Simulation of LoRaWAN LR-FHSS for Direct-to-Satellite Scenario," in IEEE Wireless Communications Letters, doi: 10.1109/LWC.2021.3135984.
%This is the function definition. It calculates the time a packet spends in the air when being transmitted.
%nputs:
%Payload: The size of the payload in bytes (the actual data being sent).
%Header_ToA_DR8: The total time for the header part of the message.
%M: The mode (2 for DR8, 4 for DR9). This determines how the payload is divided into fragments.
function [Payload_CRC_ToA_DR8,Payload_CRC_ToA_DR8_WH] = ToA_Packets_DR8(Payload,Header_ToA_DR8,M)    
   %Outputs:
%Payload_CRC_ToA_DR8: ToA for the entire packet (header + payload).
%Payload_CRC_ToA_DR8_WH: ToA for just the payload.

%Loop Over Payload Sizes
    for PL=1:length(Payload)%Loops through each payload size in the input Payload array
    Payload_CRC_ToA_DR8(PL) = Header_ToA_DR8  + ceil((Payload(PL) + 2)/M)*(102/1000); 
    %Calculates the total Time on Air for the current payload size, including the header.
    Payload_CRC_ToA_DR8_WH(PL) = ceil((Payload(PL) + 2)/M)*(102/1000); 
    end
end
