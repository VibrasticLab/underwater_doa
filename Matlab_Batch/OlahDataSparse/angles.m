function [angles]=angles(delay,c,d,d2)
% c     = kecepatan suara
% delay = waktu delay 2 mikrofon
% d     = jarak 2 mikrofon
if delay < 0
angle1=real((acosd(c*abs(delay)/d)));
angles=atand(tan(angle1*pi/180)*d/d2);
angles=180-angles;
else delay > 0;
angle1=real((acosd(c*abs(delay)/d)));
angles=atand(tan(angle1*pi/180)*d/d2);
angles = 0+angles;
end
end
