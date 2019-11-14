function [Points, G, min_xs, min_ys] = drawVolcano(r_vent, height_vent)

    H = height_vent;
    H2 = height_vent-200;
    R = 200;
    alpha = tan(0.5);
    gamma = deg2rad(50);
    w = R/(cos(alpha));
    Pa = H*tan(gamma);
    Pb = Pa*sin(alpha);
    Pc = Pa*cos(alpha);
    
    G = [0.55 0.3 0];

    P1x = R;             P1y = -R*tan(alpha);           P1z = H2;
    P2x = R;             P2y = R*tan(alpha);            P2z = H2;
    P3x = R*tan(alpha);  P3y = R;                       P3z = H2;
    P4x = -R*tan(alpha); P4y = R;                       P4z = H;
    P5x = -R;            P5y = R*tan(alpha);            P5z = H;
    P6x = -R;            P6y = -R*tan(alpha);           P6z = H2;
    P7x = -R*tan(alpha); P7y = -R;                      P7z = H2;
    P8x = R*tan(alpha);  P8y = -R;                      P8z = H2;

    P1 = [P1x P1y P1z];
    P2 = [P2x P2y P2z];
    P3 = [P3x P3y P3z];
    P4 = [P4x P4y P4z];
    P5 = [P5x P5y P5z];
    P6 = [P6x P6y P6z];
    P7 = [P7x P7y P7z];
    P8 = [P8x P8y P8z];

    P9 = [P1x+Pc P1y-Pb 0];
    P10 = [P2x+Pc P2y+Pb 0];
    P11 = [P3x+Pb P3y+Pc 0];
    P12 = [P4x-Pb P4y+Pc 0];
    P13 = [P5x-Pc P5y+Pb 0];
    P14 = [P6x-Pc P6y-Pb 0];
    P15 = [P7x-Pb P7y-Pc 0];
    P16 = [P8x+Pb P8y-Pc 0];
    
    Points.P1 = P1;    Points.P5 = P5;
    Points.P2 = P2;    Points.P6 = P6;
    Points.P3 = P3;    Points.P7 = P7;
    Points.P4 = P4;    Points.P8 = P8;
    
    Points.P9 = P9;    Points.P13 = P13;
    Points.P10 = P10;  Points.P14 = P14;
    Points.P11 = P11;  Points.P15 = P15;
    Points.P12 = P12;  Points.P16 = P16;
    
    xs = [P9(1) P10(1) P11(1) P12(1) P13(1) P14(1) P15(1) P16(1)];
    ys = [P9(2) P10(2) P11(2) P12(2) P13(2) P14(2) P15(2) P16(2)];
    
    min_xs = min(xs);
    min_ys = min(ys);
    
%     fill3([P8(1) P1(1) P9(1) P16(1)],[P8(2) P1(2) P9(2) P16(2)],[P8(3) P1(3) P9(3) P16(3)],1);
%     fill3([P1(1) P2(1) P10(1) P9(1)],[P1(2) P2(2) P10(2) P9(2)],[P1(3) P2(3) P10(3) P9(3)],1);
%     fill3([P2(1) P3(1) P11(1) P10(1)],[P2(2) P3(2) P11(2) P10(2)],[P2(3) P3(3) P11(3) P10(3)],1);
%     fill3([P3(1) P4(1) P12(1) P11(1)],[P3(2) P4(2) P12(2) P11(2)],[P3(3) P4(3) P12(3) P11(3)],1);
%     fill3([P4(1) P5(1) P13(1) P12(1)],[P4(2) P5(2) P13(2) P12(2)],[P4(3) P5(3) P13(3) P12(3)],1);
%     fill3([P5(1) P6(1) P14(1) P13(1)],[P5(2) P6(2) P14(2) P13(2)],[P5(3) P6(3) P14(3) P13(3)],1);
%     fill3([P6(1) P7(1) P15(1) P14(1)],[P6(2) P7(2) P15(2) P14(2)],[P6(3) P7(3) P15(3) P14(3)],1);
%     fill3([P7(1) P8(1) P16(1) P15(1)],[P7(2) P8(2) P16(2) P15(2)],[P7(3) P8(3) P16(3) P15(3)],1);

end