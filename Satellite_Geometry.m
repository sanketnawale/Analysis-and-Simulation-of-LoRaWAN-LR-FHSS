% Paper title: Analysis and Simulation of LoRaWAN LR-FHSS for Direct-to-Satellite Scenario
% IEEE XPlore: https://ieeexplore.ieee.org/document/9653679
% Authors: Muhammad Asad Ullah, Konstantin Mikhaylov, Hirley Alves

% Cite this: M. A. Ullah, K. Mikhaylov and H. Alves, "Analysis and Simulation of LoRaWAN LR-FHSS for Direct-to-Satellite Scenario," in IEEE Wireless Communications Letters, doi: 10.1109/LWC.2021.3135984.
%This is the function definition. It calculates key distances and angles for a satellite:
%Distance: How far the user (on the ground) is from the satellite.
%Elevation_Angles: How high the satellite appears in the sky.
%%Ground_distance: The horizontal distance on the ground under the satellite's coverage area.
%FootPrint_R: The satellite’s total coverage radius on the ground.
%inputs:

%H: The height of the satellite above the Earth.
%E: Elevation angles (the angle between the satellite and the horizon from the user’s view).
function [Distance, Elevation_Angles, Ground_distance,FootPrint_R] = Satellite_Geometry (H,E)

     %% Maximum distance from user to satellite 
    % N. Okati and T. Riihonen, "Stochastic Analysis of Satellite Broadband by Mega-Constellations with Inclined LEOs," ...
    .... 2020 IEEE 31st Annual International Symposium on Personal, Indoor and Mobile Radio Communications, London, United Kingdom, 2020, pp. 1-6
 
    % https://ieeexplore.ieee.org/document/9217379

    R = 6378e3;            %Radius of earth 

    X = cosd(E).*cosd(E);  %cos^2 %It calculates cos2(Elevation Angle)cos 2 (Elevation Angle) for all elevation angles (E).
    V = ((H + R)./R)^2;%Define a RatioThis ratio helps calculate how the satellite’s height affects the geometry of the distances.

    Slant_Range = R.*(sqrt(V-X) - sind(E));  % Slant_Range(E): distance from user to satellite
    %This calculates the slant range, which is the straight-line distance between the user and the satellite.
    %Sorts the slant range distances from smallest to largest.
    %Sorting distance 780 to 2325 km
    Sort_SR=sort( Slant_Range);  
 
    %% Points (uniform step sizes) for base of coverage triangle 
    
    %b = sqrt(c^2 - a^2)  Right angle triangle formula

    % dMAX is the maximum propogation link (user to satellite)

    FootPrint_R = sqrt(max( Slant_Range)^2 - H^2);%Calculates the satellite’s coverage radius on the ground (how far its signal reaches).
    
    Ground_distance = 1:90e3:FootPrint_R; %Step size = 90e3 km
    %Creates a list of distances on the ground from the satellite’s center, with steps of 90 km.
%Starts at 1 km and goes up to the satellite’s footprint radius

    Distance = zeros(1,length(Ground_distance));
    %Calculates the straight-line distance from the user to the satellite for each ground distance.
    for loop=1:length(Ground_distance)
    
        Distance(1,loop) = sqrt(H^2 + Ground_distance(loop)^2);
    
    end
    
    %% Finding Elevation angle as function of maximum distance
    %https://www.researchgate.net/publication/220099290 ...
    ..._The_Range_and_Horizon_Plane_Simulation_for_Ground_Stations_of_Low_Earth_Orbiting_LEO_Satellites
    
    %Equation 8
    %Calculate Elevation Angles Calculates the elevation angle for each distance using satellite geo
    Elevation_Angle_Ground_Distance=zeros(1,length(Distance));
    
    for p=1:length(Distance)
    
        
        Elevation_Angle_Ground_Distance(1,p) = (H*(H+2*R) - Distance(p)^2)/(2*Distance(p)*R);
    
    end

    Elevation_Angles=asind(Elevation_Angle_Ground_Distance);%Converts the result to degrees for easier interpretatio

end
