clc,clear,close all
% THE NUMBER OF FINITE ELEMENTS IN EACH AXIS FOR MODEL DOMAIN
Domain_nelx = 400; Domain_nely = 200;
% THE HEIGHT OF DESIGN AND SUBSTRATE DOMAIN
Design_Height = 15; Substrate_Height = 20;
% THE TARGET POINT AND INCIDENT WAVELENGTH
focus_point = [200,120]; lambda = 35;
% THE ELEMENTS INDEX OF DESIGN AND SUBSTRATE DOMAIN
Design_Element_Index = repmat(1:Domain_nely:Domain_nelx*Domain_nely,Design_Height,1);
Design_Element_Index = Design_Element_Index + ....
    repmat((Domain_nely - (Design_Height + Substrate_Height):Domain_nely - Substrate_Height - 1)',1,Domain_nelx);
Substrate_Element_Index = repmat(1:Domain_nely:Domain_nelx*Domain_nely,Substrate_Height,1);
Substrate_Element_Index = Substrate_Element_Index + ....
    repmat((Domain_nely - Substrate_Height:Domain_nely-1)',1,Domain_nelx);
% MAIN FUNCTION CALL
Top_MFSE_EM(focus_point,Design_Element_Index,Substrate_Element_Index,Domain_nelx,Domain_nely,1,3,lambda);