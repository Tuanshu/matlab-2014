clear all
Axial_resolution=9;

r=0:0.1:5500;    %micron
R=9000;      %micron
conic=0;
Cornea_function=R-((r.^2)/R)./(1+(1-(1+conic).*(r/R).^2).^0.5);

for j=1:length(r)
    if isempty(find(Cornea_function<(Cornea_function(j)-Axial_resolution),1,'first'))
        Lateral_resolution(j)=Lateral_resolution(j-1);
    else
        Lateral_resolution(j)=r(find(Cornea_function<(Cornea_function(j)-Axial_resolution),1,'first'))-r(j);
    end
end

%%
N=1;
p=1;
Cornea_function_2=R-Cornea_function;
Cornea_function_2_pre=0;
while(isempty(find(Cornea_function_2>(Cornea_function_2_pre+N),1,'first'))==0)
    index_array(p)=find(Cornea_function_2>(Cornea_function_2_pre+N),1,'first');
    Cornea_function_2_pre=Cornea_function_2(index_array(p));
    p=p+1;
end
r_array=r(index_array);
plot(r_array, '.');
xlabel('Frame Index');
ylabel('Radius of Ring (micron)');
%%
plot(r,Lateral_resolution);
xlabel('Radial Position (micron)');
%ylabel('Fringe Spacing (micron)');
ylabel('Width of Ring (micron)');