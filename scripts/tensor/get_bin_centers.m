function [x,y,z,Vx,Vy,Vz] = get_bin_centers( json );
x  = [json.minval(1): json.binwidth(1) : json.maxval(1)];
y  = [json.minval(2): json.binwidth(2) : json.maxval(2)];
z  = [json.minval(3): json.binwidth(3) : json.maxval(3)];
Vx = [json.minval(4): json.binwidth(4) : json.maxval(4)];
Vy = [json.minval(5): json.binwidth(5) : json.maxval(5)];
Vz = [json.minval(6): json.binwidth(6) : json.maxval(6)];
