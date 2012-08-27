for z = 1:62
    filename = sprintf('piece%02d.txt',z);
    Bez_Spline(filename,z)
end

%d = dir('piece*.txt');
%for z = 10:62
    %filename = d(z).name
    %Bez_Spline(filename,z)
%end

