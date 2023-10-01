function y=f1(f,dmu)
if abs(f+dmu)<=1
    y=f+dmu;
else
    y=sign(f+dmu);
end
end