function str = TFCFormat(s)

     str = [];

for j=1:6
     [t,s] = strtok(s,',');
     str=[str t ','];
end

for j=1:2

[t,s] = strtok(s,',');

if(isempty(strfind(t,'.')))
     str=[str t '.0,'];
   else
     str=[str t ','];
end

end

[t,s] = strtok(s,',');

if(isempty(strfind(t,'.')))
     str=[str t '.0'];
   else
     str=[str t];
end

str = [str s];

